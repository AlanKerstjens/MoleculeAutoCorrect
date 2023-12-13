#ifndef _MOLECULE_AUTO_CORRECT_HPP_
#define _MOLECULE_AUTO_CORRECT_HPP_

#include "TreeSearch.hpp"
#include "MoleculePerturber.hpp"
#include "MolecularPerturbationUtils.hpp"
#include "FunctionalArithmetic.hpp"
#include "ChemicalDictionary.hpp"

namespace MoleculeAutoCorrect {

MoleculePerturber DefaultMoleculePerturber() {
  MoleculePerturber perturber;
  // Because the vertices of the tree search are RDKit::ROMol, meaning that they
  // aren't modified, the SSSR and topological distance matrix aren't invalidated,
  // and it's worthwhile pre-computing them and using them during perturbations.
  perturber.assess_connectivity_with_sssr = true;
  perturber.assess_distances_with_distance_matrix = true;
  return perturber;
};

struct Settings {
  const MoleculePerturber* perturber = nullptr;
  const ChemicalDictionary* dictionary = nullptr;
  const MolecularConstraints* constraints = nullptr;

  // By default the least significant perturbations, that is, "decorations", are
  // applied before the more significant perturbations changing the topology.
  // Deletions are prioritized over insertions as they simplify molecules,
  // making it more likely that they lead to a correct molecule.
  std::vector<MolecularPerturbation::Type> perturbation_type_priority {
    MolecularPerturbation::Type::ExplicitHydrogensChange_t,
    MolecularPerturbation::Type::FormalChargeChange_t,
    MolecularPerturbation::Type::AtomicNumberChange_t,
    MolecularPerturbation::Type::BondTypeChange_t,
    MolecularPerturbation::Type::BondDeletion_t,
    MolecularPerturbation::Type::AtomDeletion_t,
    MolecularPerturbation::Type::BondInsertion_t,
    MolecularPerturbation::Type::AtomInsertion_t
  };

  bool sanitize_after_expansion = true;

  bool allow_degree_change_during_valence_correction = true;
  // Set settings below to false if you are experiencing alkane growth.
  bool attempt_degree_correction_with_insertions = false;
  bool attempt_bond_correction_with_atom_insertions = true;
  bool attempt_environment_correction_with_atom_insertions = true;

  bool constrain_z_on_dv = true;
  bool constrain_q_on_dvz = true;
  bool constrain_h_on_dvzq = true;
  bool constrain_b_on_k1k2 = true;

  std::size_t max_tree_size = 25000;
  std::size_t max_tree_depth = 25;
  std::size_t n_solutions = 1;
  std::size_t n_top_solutions = 1;
};

std::vector<CircularAtomicEnvironment> ForeignEnvironments(
  const RDKit::ROMol& molecule,
  const ChemicalDictionary& dictionary) {
  std::vector<CircularAtomicEnvironment> foreign_environments;
  foreign_environments.reserve(molecule.getNumAtoms());
  for (const RDKit::Atom* atom : molecule.atoms()) {
    CircularAtomicEnvironment environment = dictionary.Environment(atom);
    if (dictionary.IsForeignEnvironment(environment.Key())) {
      foreign_environments.push_back(std::move(environment));
    };
  };
  return foreign_environments;
};

std::vector<MolecularPerturbation::Target> TargetPriority(
  const RDKit::ROMol& molecule,
  const std::vector<AtomKey>& atom_keys,
  const std::vector<BondKey>& bond_keys,
  const std::vector<CircularAtomicEnvironment>& foreign_environments,
  const ChemicalDictionary& dictionary) {
  std::vector<MolecularPerturbation::Target> priorities;
  if (foreign_environments.empty()) {
    return priorities;
  };
  auto [atom_overlaps, bond_overlaps] =
    CircularEnvironmentOverlap(foreign_environments);
  std::size_t n_atoms = atom_keys.size(), n_bonds = bond_keys.size();
  priorities.reserve(n_atoms + n_bonds);
  for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
    // We don't care about atoms or bonds with null overlap.
    if (atom_overlaps[atom_idx] > 0) {
      priorities.emplace_back(MolecularPerturbation::TargetType::Atom, atom_idx);
    };
  };
  for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
    if (bond_overlaps[bond_idx] > 0) {
      priorities.emplace_back(MolecularPerturbation::TargetType::Bond, bond_idx);
    };
  };
  // Sorts graph elements (atoms and bonds) by their overlap in descending order.
  // Tied elements are sorted by their normalized frequencies in the Chemical
  // Dictionary in ascending order.
  // Future improvement: Currently there is a relatively high chance of ties. 
  // One could reduce the probability of ties by calculating the overlap of 
  // arbitrary subgraphs as opposed to circular environments.
  double total_atom_frequency = dictionary.GetTotalAtomFrequency();
  double total_bond_frequency = dictionary.GetTotalBondFrequency();
  std::sort(priorities.begin(), priorities.end(),
    [&] (const auto& p1, const auto& p2) {
      unsigned o1 = p1.first == MolecularPerturbation::TargetType::Atom ?
        atom_overlaps[p1.second] :
        bond_overlaps[p1.second];
      unsigned o2 = p2.first == MolecularPerturbation::TargetType::Atom ?
        atom_overlaps[p2.second] :
        bond_overlaps[p2.second];
      if (o1 != o2) {
        return o1 > o2;
      };
      double nf1, nf2;
      if (p1.first == MolecularPerturbation::TargetType::Atom) {
        nf1 = dictionary.AtomFrequency(
          atom_keys[p1.second]) / total_atom_frequency;
      } else { // p1.first == MolecularPerturbation::TargetType::Bond
        // Given that BondKeys are identified by a pair of AtomKeys, they are
        // intrinsically more unique and less frequent. Assuming a random
        // distribution, if the probability of encountering a duplicate AtomKey
        // is P, that of encountering a duplicate BondKey approaches P^2.
        double bond_frequency = dictionary.BondFrequency(bond_keys[p1.second]);
        nf1 = bond_frequency * bond_frequency / total_bond_frequency;
      };
      if (p2.first == MolecularPerturbation::TargetType::Atom) {
        nf2 = dictionary.AtomFrequency(
          atom_keys[p2.second]) / total_atom_frequency;
      } else {
        double bond_frequency = dictionary.BondFrequency(bond_keys[p2.second]);
        nf2 = bond_frequency * bond_frequency / total_bond_frequency;
      };
      return nf1 < nf2;
    }
  );
  return priorities;
};

class Vertex : public RDKit::ROMol {
  const Settings* settings = nullptr;
  MolecularPerturbationQueue queue;
  MolecularKeys molecular_keys;
  std::vector<CircularAtomicEnvironment> foreign_environments;
  std::vector<MolecularPerturbation::Target> perturbation_targets;
  std::vector<MolecularPerturbation::Type> perturbation_types;
  std::size_t perturbation_target_idx = 0, perturbation_type_idx = 0;
  std::size_t n_foreign_atoms = 0;
  std::size_t n_foreign_bonds = 0;
  std::size_t n_foreign_environments = 0;
  bool enqueue_perturbations = false;
  bool enqueued_atom_key_correcting_perturbations = false;
  bool enqueued_bond_key_correcting_perturbations = false;
  bool setup_environment_perturbations = false;

public:
  Vertex() = default; // Default constructor to store Vertices in vectors.
  Vertex(
    const RDKit::ROMol& molecule,
    const Settings* settings) :
    RDKit::ROMol(molecule),
    settings(settings),
    molecular_keys(molecule) {
    TagMolecule(*this, true);
    foreign_environments = ForeignEnvironments(*this, *settings->dictionary);
    n_foreign_atoms = NForeignAtoms();
    n_foreign_bonds= NForeignBonds();
    n_foreign_environments = foreign_environments.size();
    if (!IsFamiliar()) {
      enqueue_perturbations = true;
    };
  };

private:
  void EnqueueDegreeCorrectingPerturbatons(std::size_t atom_idx) {
    settings->perturber->AtomDeletions(
      queue, *this, atom_idx, settings->constraints);
    for (const RDKit::Bond* bond : atomBonds(getAtomWithIdx(atom_idx))) {
      // Deletions of neighboring atoms (with reconnection of remaining 
      // neighbors) can only reduce the degree if the neighbor is a peripheral 
      // atom. We allow the deletion of non-peripheral neighbors anyways to 
      // allow for the deletion of entire branches of the molecule over the 
      // span of multiple expansions.
      settings->perturber->AtomDeletions(
        queue, *this, bond->getOtherAtomIdx(atom_idx), settings->constraints);
      settings->perturber->BondDeletions(
        queue, *this, bond->getIdx(), settings->constraints);
    };
    if (settings->attempt_degree_correction_with_insertions) {
      settings->perturber->AtomInsertions(
        queue, *this, atom_idx, settings->constraints);
      settings->perturber->BondInsertions(
        queue, *this, atom_idx, settings->constraints);
    };
  };

  void EnqueueValenceCorrectingPerturbations(std::size_t atom_idx) {
    for (const RDKit::Bond* bond : atomBonds(getAtomWithIdx(atom_idx))) {
      settings->perturber->BondTypeChanges(
        queue, *this, bond->getIdx(), settings->constraints);
    };
    if (settings->allow_degree_change_during_valence_correction) {
      EnqueueDegreeCorrectingPerturbatons(atom_idx);
    };
  };

  void EnqueueAtomicNumberChanges(std::size_t atom_idx) {
    if (settings->constrain_z_on_dv) {
      // It is tempting to check if the current value is already allowed and, 
      // if that is the case, skip enqueing perturbations. However, we can't 
      // guarantee that it is possible to fix the key/environment by keeping it.
      // Conversely, there may be other valid values that allow us to fix the 
      // key/environment.
      const AtomKey& atom_key = molecular_keys.atom_keys[atom_idx];
      std::vector<std::uint8_t> allowed_atomic_numbers = 
        settings->dictionary->Z_DV(atom_key.dv());
      settings->perturber->AtomicNumberChanges(
        queue, *this, atom_idx, settings->constraints, &allowed_atomic_numbers);
    } else {
      settings->perturber->AtomicNumberChanges(
        queue, *this, atom_idx, settings->constraints);
    };
  };

  void EnqueueFormalChargeChanges(std::size_t atom_idx) {
    if (settings->constrain_q_on_dvz) {
      const AtomKey& atom_key = molecular_keys.atom_keys[atom_idx];
      std::vector<std::int8_t> allowed_formal_charges = 
        settings->dictionary->Q_DVZ(atom_key.dvz());
      settings->perturber->FormalChargeChanges(
        queue, *this, atom_idx, settings->constraints, &allowed_formal_charges);
    } else {
      settings->perturber->FormalChargeChanges(
        queue, *this, atom_idx, settings->constraints);
    };
  };

  void EnqueueExplicitHydrogenChanges(std::size_t atom_idx) {
    if (settings->constrain_h_on_dvzq) {
      const AtomKey& atom_key = molecular_keys.atom_keys[atom_idx];
      std::vector<std::uint8_t> allowed_n_explicit_hydrogens = 
        settings->dictionary->H_DVZQ(atom_key.dvzq());
      settings->perturber->ExplicitHydrogenChanges(
        queue, *this, atom_idx,
        settings->constraints, &allowed_n_explicit_hydrogens);
    } else {
      settings->perturber->ExplicitHydrogenChanges(
        queue, *this, atom_idx, settings->constraints);
    };
  };

  void EnqueueBondTypeChanges(std::size_t bond_idx) {
    if (settings->constrain_b_on_k1k2) {
      const BondKey& bond_key = molecular_keys.bond_keys[bond_idx];
      std::vector<RDKit::Bond::BondType> allowed_bond_types = 
        settings->dictionary->B_K1K2(bond_key.k1k2());
      settings->perturber->BondTypeChanges(
        queue, *this, bond_idx, settings->constraints, &allowed_bond_types);
    } else {
      settings->perturber->BondTypeChanges(
        queue, *this, bond_idx, settings->constraints);     
    };
  };

  void EnqueueAtomKeyCorrectingPerturbations() {
    std::size_t n_atoms = molecular_keys.atom_keys.size();
    for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
      const AtomKey& atom_key = molecular_keys.atom_keys[atom_idx];
      // AtomKeys and their partial keys are considered foreign when their
      // frequency in the ChemicalDictionary is lower than the Settings threshold.
      AtomKey::Error error = settings->dictionary->AtomKeyError(atom_key);
      switch (error) {
        case AtomKey::Error::D:
          // If a partial key is foreign all higher order partial keys are also
          // foreign. There's no point in trying to fix higher order properties
          // before the lower order one has been corrected.
          EnqueueDegreeCorrectingPerturbatons(atom_idx);
          continue;
        case AtomKey::Error::DV:
          EnqueueValenceCorrectingPerturbations(atom_idx);
          continue;
        case AtomKey::Error::DVZ:
          EnqueueAtomicNumberChanges(atom_idx);
          continue;
        case AtomKey::Error::DVZQ:
          EnqueueFormalChargeChanges(atom_idx);
          continue;
        case AtomKey::Error::DVZQH:
          EnqueueExplicitHydrogenChanges(atom_idx);
      };
    };
  };

  void EnqueueBondKeyCorrectingPerturbations() {
    std::size_t n_bonds = molecular_keys.bond_keys.size();
    for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
      const BondKey& bond_key = molecular_keys.bond_keys[bond_idx];
      const RDKit::Bond* bond = getBondWithIdx(bond_idx);
      std::size_t begin_atom_idx = bond->getBeginAtomIdx();
      std::size_t end_atom_idx = bond->getEndAtomIdx();
      BondKey::Error error = settings->dictionary->BondKeyError(bond_key);
      switch (error) {
        case BondKey::Error::K1K2:
          EnqueueExplicitHydrogenChanges(begin_atom_idx);
          EnqueueExplicitHydrogenChanges(end_atom_idx);
          EnqueueFormalChargeChanges(begin_atom_idx);
          EnqueueFormalChargeChanges(end_atom_idx);
          EnqueueAtomicNumberChanges(begin_atom_idx);
          EnqueueAtomicNumberChanges(end_atom_idx);
          settings->perturber->BondDeletions(
            queue, *this, bond_idx, settings->constraints);
          settings->perturber->AtomDeletions(
            queue, *this, begin_atom_idx, settings->constraints);
          settings->perturber->AtomDeletions(
            queue, *this, end_atom_idx, settings->constraints);
          if (settings->attempt_bond_correction_with_atom_insertions) {
            settings->perturber->AtomInsertions(
              queue, *this, begin_atom_idx, settings->constraints);
            settings->perturber->AtomInsertions(
              queue, *this, end_atom_idx, settings->constraints);
          };
          settings->perturber->BondInsertions(
            queue, *this, begin_atom_idx, settings->constraints);
          settings->perturber->BondInsertions(
            queue, *this, end_atom_idx, settings->constraints);
          continue;
        case BondKey::Error::K1K2B:
          EnqueueBondTypeChanges(bond_idx);
      };
    };
  };

  void AdvanceEnvironmentPerturbationIndices() {
    if (++perturbation_type_idx >= perturbation_types.size()) {
      perturbation_type_idx = 0;
      if (++perturbation_target_idx >= perturbation_targets.size()) {
        enqueue_perturbations = false;
      };
    };
  };

  void EnqueueEnvironmentCorrectingPerturbations() {
    MolecularPerturbation::Type perturbation_type =
      perturbation_types[perturbation_type_idx];
    auto [target_type, target_idx] =
      perturbation_targets[perturbation_target_idx];
    if (target_type == MolecularPerturbation::TargetType::Atom) {
      switch (perturbation_type) {
        case MolecularPerturbation::Type::ExplicitHydrogensChange_t:
          EnqueueExplicitHydrogenChanges(target_idx);
          break;
        case MolecularPerturbation::Type::FormalChargeChange_t:
          EnqueueFormalChargeChanges(target_idx);
          break;
        case MolecularPerturbation::Type::AtomicNumberChange_t:
          EnqueueAtomicNumberChanges(target_idx);
          break;
        case MolecularPerturbation::Type::AtomDeletion_t:
          settings->perturber->AtomDeletions(
            queue, *this, target_idx, settings->constraints);
          break;
        case MolecularPerturbation::Type::AtomInsertion_t:
          if (settings->attempt_environment_correction_with_atom_insertions) {
            settings->perturber->AtomInsertions(
              queue, *this, target_idx, settings->constraints);
          };
          break;
        case MolecularPerturbation::Type::BondInsertion_t:
          settings->perturber->BondInsertions(
            queue, *this, target_idx, settings->constraints);
      };
    } else if (target_type == MolecularPerturbation::TargetType::Bond) {
      switch (perturbation_type) {
        case MolecularPerturbation::Type::BondTypeChange_t:
          EnqueueBondTypeChanges(target_idx);
          break;
        case MolecularPerturbation::Type::BondDeletion_t:
          settings->perturber->BondDeletions(
            queue, *this, target_idx, settings->constraints);
      };
    };
    AdvanceEnvironmentPerturbationIndices();
  };

  void EnqueuePerturbations() {
    if (!enqueue_perturbations) {
      return;
    };
    if (n_foreign_atoms && !enqueued_atom_key_correcting_perturbations) {
      EnqueueAtomKeyCorrectingPerturbations();
      enqueued_atom_key_correcting_perturbations = true;
      // Since BondKeys depend on AtomKeys, if there are foreign AtomKeys the
      // corresponding BondKeys are also foreign. It's pointless to try to fix
      // these BondKeys, as the foreign AtomKeys are a more pressing issue.
      // Note however that the absence of foreign AtomKeys does NOT guarantee
      // the absence of foreign BondKeys.
      return;
    };
    if (n_foreign_bonds && !enqueued_bond_key_correcting_perturbations) {
      EnqueueBondKeyCorrectingPerturbations();
      enqueued_bond_key_correcting_perturbations = true;
      // Since environments depend on keys, if there are foreign keys the
      // corresponding environments are also foreign. The absence of foreign
      // keys doesn't guarantee the absence of foreign environments.
      return;
    };
    if (n_foreign_environments) {
      if (!setup_environment_perturbations) {
        const MolecularPerturbation::TypeMask& ptm =
          settings->perturber->perturbation_types;
        perturbation_types.reserve(ptm.count());
        for (MolecularPerturbation::Type type :
          settings->perturbation_type_priority) {
          if (ptm[type]) {
            perturbation_types.push_back(type);
          };
        };
        perturbation_targets = TargetPriority(
          *this,
          molecular_keys.atom_keys,
          molecular_keys.bond_keys,
          foreign_environments,
          *(settings->dictionary));
        setup_environment_perturbations = true;
      };
      EnqueueEnvironmentCorrectingPerturbations();
      return;
    };
    enqueue_perturbations = false; // We should never get here.
  };

  std::size_t NForeignAtoms() const {
    std::size_t n_foreign_atoms = 0;
    for (const AtomKey& atom_key : molecular_keys.atom_keys) {
      if (settings->dictionary->IsForeignAtom(atom_key)) {
        ++n_foreign_atoms;
      };
    };
    return n_foreign_atoms;
  };

  std::size_t NForeignBonds() const {
    std::size_t n_foreign_bonds = 0;
    for (const BondKey& bond_key : molecular_keys.bond_keys) {
      if (settings->dictionary->IsForeignBond(bond_key)) {
        ++n_foreign_bonds;
      };
    };
    return n_foreign_bonds;
  };

public:
  std::optional<Vertex> Expand() {
    std::optional<Vertex> perturbed_vertex;
    while (queue.empty() && enqueue_perturbations) {
      EnqueuePerturbations();
    };
    if (!queue.empty()) {
      const std::shared_ptr<MolecularPerturbation>& perturbation = queue.front();
      RDKit::RWMol perturbed_molecule = (*perturbation)(*this);
      if (settings->sanitize_after_expansion) {
        PartialSanitization(perturbed_molecule);
      };
      perturbed_vertex = Vertex(std::move(perturbed_molecule), settings);
      queue.pop();
    };
    return perturbed_vertex;
  };

  bool IsExpandable() const {
    return !queue.empty() || enqueue_perturbations;
  };

  double Familiarity1() const {
    double n_foreign_keys = 
      n_foreign_atoms + n_foreign_bonds + n_foreign_environments;
    double n_keys = getNumAtoms() * 2.0 + getNumBonds();
    return (n_keys - n_foreign_keys) / n_keys;
  };

  double Familiarity2() const {
    double n_atoms = getNumAtoms();
    return (n_atoms - n_foreign_environments) / n_atoms;
  };

  bool IsFamiliar() const {
    return Familiarity1() >= 1.0;
  };
};

}; // ! MoleculeAutoCorrect namespace


// std::hash specialization for MoleculeAutoCorrect::Vertex.
template <>
class std::hash<MoleculeAutoCorrect::Vertex> {
  static inline std::hash<RDKit::ROMol> molecule_hasher;
public:
  std::size_t operator()(const MoleculeAutoCorrect::Vertex& vertex) const {
    return molecule_hasher(vertex);
  };
};


namespace MoleculeAutoCorrect {

double Familiarity(const Vertex& vertex) {
  return vertex.Familiarity1();
};

std::optional<Vertex> Expansion(Vertex& vertex) {
  return vertex.Expand();
};

struct MoleculeDiscovery {
  RDKit::ROMol molecule;
  TreeSearch<Vertex>::Tree::vertex_descriptor vertex = 0;
  std::size_t order = 0;
  std::size_t depth = 0;
  double familiarity = 0.0;
};

struct Result {
  std::size_t n_expansions = 0;
  std::vector<MoleculeDiscovery> top_discoveries;
};


namespace Policy {


namespace Objective {

// The function signature of our objectives is:
// double(
//   const TreeSearch<Vertex>&,
//   TreeSearch<Vertex>::Tree::vertex_descriptor)>
// The reason for using this signature as opposed to double(const RDKit::ROMol&)
// is that the former allows us to retrieve information about the Vertex and its
// position within the TreeSearch, enabling more sophisticated policies.

// Combining objectives is a powerful way of fine-tuning policies. The Function
// class wraps a function with the signature double(Args...) and overloads 
// arithmetic operators for it. Wrapper is a specialization of Function where
// Args... is set to match the signature of our objectives.
typedef Function<
  const TreeSearch<Vertex>&,
  TreeSearch<Vertex>::Tree::vertex_descriptor> Wrapper;

// Users might get confused by the template classes. They may prefer to define
// objectives as double(const RDKit::ROMol&). This class wraps the former as a
// double(Args...).
class MoleculeObjective {
  std::function<double(const RDKit::ROMol&)> objective;
public:
  MoleculeObjective(
    const std::function<double(const RDKit::ROMol&)>& objective) : 
      objective(objective) {};
  
  double operator()(
    const TreeSearch<Vertex>& tree_search,
    TreeSearch<Vertex>::Tree::vertex_descriptor v) const {
    return objective(tree_search.GetVertex(v));
  };
};

double Familiarity1(
  const TreeSearch<Vertex>& tree_search,
  TreeSearch<Vertex>::Tree::vertex_descriptor v) {
  return tree_search.GetVertex(v).Familiarity1();
};

double Familiarity2(
  const TreeSearch<Vertex>& tree_search,
  TreeSearch<Vertex>::Tree::vertex_descriptor v) {
  return tree_search.GetVertex(v).Familiarity2();
};

class TopologicalSimilarity {
  std::shared_ptr<RDKit::SparseIntVect<std::uint32_t>> reference_fingerprint;
  unsigned fingerprint_radius = 2;

public:
  TopologicalSimilarity(
    const RDKit::ROMol& reference_molecule,
    unsigned fingerprint_radius = 2) :
    reference_fingerprint(RDKit::MorganFingerprints::getFingerprint(
      reference_molecule, fingerprint_radius)),
    fingerprint_radius(fingerprint_radius) {};

  double operator()(
    const TreeSearch<Vertex>& tree_search,
    TreeSearch<Vertex>::Tree::vertex_descriptor v) const {
    const Vertex& vertex = tree_search.GetVertex(v);
    RDKit::SparseIntVect<std::uint32_t>* fingerprint = 
      RDKit::MorganFingerprints::getFingerprint(vertex, fingerprint_radius);
    double similarity = RDKit::TanimotoSimilarity(
      *reference_fingerprint, *fingerprint);
    delete fingerprint;
    return similarity;
  };
};

struct Constant {
  double x = 0.0;
  Constant(double x) : x(x) {};
  double operator()(
    const TreeSearch<Vertex>& tree_search,
    TreeSearch<Vertex>::Tree::vertex_descriptor v) const {
    return x;
  };
};

}; // ! MoleculeAutoCorrect::Policy::Objective namespace


// Virtually all policies can be expressed as a GreedyPolicy, provided that we
// adjust the objective on which the GreedyPolicy selects and tries to maximize.
enum class Type {
  BFS,
  Familiarity,
  DistanceNormalizedFamiliarity,
  Astar,
  UCT,
  MLR
};

struct BFS : GreedyPolicy<Vertex> {
  BFS(const RDKit::ROMol& root_molecule) : 
    GreedyPolicy<Vertex>(Objective::TopologicalSimilarity(root_molecule)) {};
};

struct Familiarity : GreedyPolicy<Vertex> {
  Familiarity() : GreedyPolicy<Vertex>(Objective::Familiarity1) {};
};

struct DistanceNormalizedFamiliarity : GreedyPolicy<Vertex> {
  DistanceNormalizedFamiliarity(const RDKit::ROMol& root_molecule) : 
    GreedyPolicy<Vertex>(
      Objective::Wrapper(Objective::Familiarity1) / 
      (1.0 - Objective::Wrapper(Objective::TopologicalSimilarity(root_molecule)))) {};
};

struct Astar : GreedyPolicy<Vertex> {
  Astar(const RDKit::ROMol& root_molecule) :
    // Selecting the vertex with the highest sum of similarities is equal to 
    // selecting the vertex with the lowest sum of distances.
    GreedyPolicy<Vertex>(
      Objective::Wrapper(Objective::TopologicalSimilarity(root_molecule)) +
      Objective::Wrapper(Objective::Familiarity1)) {};
};

typedef UpperConfidenceTree<Vertex> UCT;

struct MLR : GreedyPolicy<Vertex> {
  MLR(const RDKit::ROMol& root_molecule) : GreedyPolicy<Vertex>(
    1.0 - ( // Maximize the complement of the MLR distance (i.e. similarity)
      0.42 
      * (1.0 - Objective::Wrapper(Objective::TopologicalSimilarity(root_molecule)))
      + 0.91 
      * (1.0 - Objective::Wrapper(Objective::Familiarity1))
      + 0.27
    )) {};
};

struct ObjectivePreservation : GreedyPolicy<Vertex> {
  ObjectivePreservation(
    const std::function<double(const RDKit::ROMol&)>& objective) :
    GreedyPolicy<Vertex>(
      Objective::Wrapper(Objective::Familiarity1) * 
      Objective::Wrapper(Objective::MoleculeObjective(objective))) {};
};

struct Dummy : GreedyPolicy<Vertex> {
  Dummy(double x) : GreedyPolicy<Vertex>(Objective::Constant(x)) {};
};

}; // ! MoleculeAutoCorrect::Policy namespace


typedef TreeSearch<Vertex>::SelectionPolicy SelectionPolicy;
typedef TreeSearch<Vertex>::RewardFunction RewardFunction;

class TerminationPolicy {
  std::size_t n_familiar_vertices = 0;
  std::size_t max_n_familiar_vertices = 0;

public:
  TerminationPolicy(std::size_t max_n_familiar_vertices = 1) :
    max_n_familiar_vertices(max_n_familiar_vertices) {};

  bool operator()(const TreeSearch<Vertex>& tree_search) {
    auto lv = tree_search.GetLastVertexDescriptor();
    const Vertex& last_vertex = tree_search.GetVertex(lv);
    n_familiar_vertices += last_vertex.IsFamiliar();
    return n_familiar_vertices >= max_n_familiar_vertices;
  };
};

}; // ! MoleculeAutoCorrect namespace


MoleculeAutoCorrect::Result AutoCorrectMolecule(
  const RDKit::ROMol& molecule,
  const MoleculeAutoCorrect::Settings& settings,
  const MoleculeAutoCorrect::SelectionPolicy& selection_policy,
  const MoleculeAutoCorrect::RewardFunction& reward_function = nullptr) {
  using namespace MoleculeAutoCorrect;
  Result result;
  result.top_discoveries.reserve(settings.n_top_solutions);

  TreeSearch<Vertex> tree_search (
    Vertex(molecule, &settings),
    Expansion,
    reward_function,
    settings.max_tree_size + 1,
    settings.max_tree_depth);
  result.n_expansions = tree_search.Search(
    selection_policy,
    TerminationPolicy(settings.n_solutions)) - 1;
  
  auto top_vertices = tree_search.TopVertices(
    Familiarity, settings.n_top_solutions);
  const auto& vertex_depths = tree_search.GetVertexDepths();
  for (auto [v, familiarity] : top_vertices) {
    const Vertex& vertex = tree_search.GetVertex(v);
    result.top_discoveries.emplace_back(
      vertex, v, v, vertex_depths[v], familiarity);
  };
  return result;
};

MoleculeAutoCorrect::Result AutoCorrectMolecule(
  const RDKit::ROMol& molecule,
  const MoleculeAutoCorrect::Settings& settings,
  MoleculeAutoCorrect::Policy::Type selection_policy_type = 
    MoleculeAutoCorrect::Policy::Type::MLR) {
  using namespace MoleculeAutoCorrect;
  SelectionPolicy selection_policy = nullptr;
  RewardFunction reward_function = nullptr;
  switch (selection_policy_type) {
    case Policy::Type::BFS:
      selection_policy = Policy::BFS(molecule);
      break;
    case Policy::Type::Familiarity:
      selection_policy = Policy::Familiarity();
      break;
    case Policy::Type::DistanceNormalizedFamiliarity:
      selection_policy = Policy::DistanceNormalizedFamiliarity(molecule);
      break;
    case Policy::Type::Astar:
      selection_policy = Policy::Astar(molecule);
      break;
    case Policy::Type::UCT:
      selection_policy = Policy::UCT(0.5);
      reward_function = Familiarity;
      break;
    case Policy::Type::MLR:
      selection_policy = Policy::MLR(molecule);
      break;
    default:
      selection_policy = Policy::MLR(molecule);
  };
  return AutoCorrectMolecule(molecule, settings,
    selection_policy, reward_function);
};

#endif // !_MOLECULE_AUTO_CORRECT_HPP_

#ifndef _PY_MOLECULE_AUTO_CORRECT_HPP_
#define _PY_MOLECULE_AUTO_CORRECT_HPP_

#include "MoleculeAutoCorrect.hpp"
#include "pySTL.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

python::list GetPerturbationTypePriority(
  const MoleculeAutoCorrect::Settings& settings) {
  return to_list(settings.perturbation_type_priority);
};

void SetPerturbationTypePriority(
  MoleculeAutoCorrect::Settings& settings,
  const python::list& new_perturbation_type_priority) {
  settings.perturbation_type_priority = 
    to_vector<MolecularPerturbation::Type>(new_perturbation_type_priority);
};

python::list GetTopDiscoveries(
  const MoleculeAutoCorrect::Result& result) {
  return to_list(result.top_discoveries);
};

void WrapMoleculeAutoCorrect() {
  python::scope module_scope = python::scope();

  boost::python::pointer_wrapper<const MolecularConstraints*>
    null_constraints (nullptr);

  python::def(
    "DefaultMoleculePerturber", MoleculeAutoCorrect::DefaultMoleculePerturber);

  python::scope settings_scope = 
  python::class_<MoleculeAutoCorrect::Settings>("Settings", 
    python::init<
      const MoleculePerturber*, 
      const ChemicalDictionary*, 
      const MolecularConstraints*>((
      python::arg("perturber"),
      python::arg("dictionary"),
      python::arg("constraints") = null_constraints)))
    
    .def_readwrite("perturber",
      &MoleculeAutoCorrect::Settings::perturber)
    .def_readwrite("dictionary",
      &MoleculeAutoCorrect::Settings::dictionary)
    .def_readwrite("constraints",
      &MoleculeAutoCorrect::Settings::constraints)
    .def_readwrite("tree_policy_type",
      &MoleculeAutoCorrect::Settings::tree_policy_type)
    .def_readwrite("uct_c",
      &MoleculeAutoCorrect::Settings::uct_c)
    .add_property("perturbation_type_priority", 
      GetPerturbationTypePriority, SetPerturbationTypePriority)
    .def_readwrite("sanitize_after_expansion",
      &MoleculeAutoCorrect::Settings::sanitize_after_expansion)
    .def_readwrite("allow_degree_change_during_valence_correction", 
      &MoleculeAutoCorrect::Settings::allow_degree_change_during_valence_correction)
    .def_readwrite("attempt_degree_correction_with_insertions", 
      &MoleculeAutoCorrect::Settings::attempt_degree_correction_with_insertions)
    .def_readwrite("attempt_bond_correction_with_atom_insertions", 
      &MoleculeAutoCorrect::Settings::attempt_bond_correction_with_atom_insertions)
    .def_readwrite("attempt_environment_correction_with_atom_insertions", 
      &MoleculeAutoCorrect::Settings::attempt_environment_correction_with_atom_insertions)
    .def_readwrite("constrain_z_on_dv", 
      &MoleculeAutoCorrect::Settings::constrain_z_on_dv)
    .def_readwrite("constrain_q_on_dvz", 
      &MoleculeAutoCorrect::Settings::constrain_q_on_dvz)
    .def_readwrite("constrain_h_on_dvzq", 
      &MoleculeAutoCorrect::Settings::constrain_h_on_dvzq)
    .def_readwrite("constrain_b_on_k1k2", 
      &MoleculeAutoCorrect::Settings::constrain_b_on_k1k2)
    .def_readwrite("max_tree_size", 
      &MoleculeAutoCorrect::Settings::max_tree_size)
    .def_readwrite("max_tree_depth", 
      &MoleculeAutoCorrect::Settings::max_tree_depth);
  
  python::enum_<MoleculeAutoCorrect::Settings::TreePolicyType>("TreePolicyType")
    .value("Familiarity", 
      MoleculeAutoCorrect::Settings::TreePolicyType::Familiarity)
    .value("BFS", 
      MoleculeAutoCorrect::Settings::TreePolicyType::BFS)
    .value("DistanceNormalizedFamiliarity", 
      MoleculeAutoCorrect::Settings::TreePolicyType::DistanceNormalizedFamiliarity)
    .value("UCT", 
      MoleculeAutoCorrect::Settings::TreePolicyType::UCT)
    .value("Astar", 
      MoleculeAutoCorrect::Settings::TreePolicyType::Astar)
    .value("MLR", 
      MoleculeAutoCorrect::Settings::TreePolicyType::MLR);
  
  python::scope mod_scope (module_scope);

  python::class_<MoleculeAutoCorrect::MoleculeDiscovery>(
    "MoleculeDiscovery", python::no_init)
    .def_readonly("molecule",
      &MoleculeAutoCorrect::MoleculeDiscovery::molecule)
    .def_readonly("order",
      &MoleculeAutoCorrect::MoleculeDiscovery::order)
    .def_readonly("depth",
      &MoleculeAutoCorrect::MoleculeDiscovery::depth)
    .def_readonly("familiarity",
      &MoleculeAutoCorrect::MoleculeDiscovery::familiarity);

  python::class_<MoleculeAutoCorrect::Result>(
    "Result", python::no_init)
    .def_readonly("n_expansions",
      &MoleculeAutoCorrect::Result::n_expansions)
    .add_property("top_discoveries", 
      GetTopDiscoveries);

  python::def("AutoCorrectMolecule", AutoCorrectMolecule, (
    python::arg("molecule"),
    python::arg("settings"),
    python::arg("n_solutions") = 1,
    python::arg("n_top_solutions") = 1));
};

#endif //_PY_MOLECULE_AUTO_CORRECT_HPP_
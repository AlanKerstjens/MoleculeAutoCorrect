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

class PythonMoleculeObjective {
  python::object objective;

public:
  PythonMoleculeObjective(const python::object& objective) : 
    objective(objective) {
    if (objective.is_none() || !PyCallable_Check(objective.ptr())) {
      PyErr_SetString(PyExc_TypeError, "Expected a callable object");
      throw boost::python::error_already_set();
    };
  };

  double operator()(const RDKit::ROMol& molecule) const {
    return python::call<double>(objective.ptr(), boost::cref(molecule));
  };
};

MoleculeAutoCorrect::Result AutoCorrectMolecule(
  const RDKit::ROMol& molecule,
  const MoleculeAutoCorrect::Settings& settings,
  const python::object& objective) {
  return AutoCorrectMolecule(molecule, settings,
    MoleculeAutoCorrect::Policy::ObjectivePreservation(
      PythonMoleculeObjective(objective)));
};


void WrapMoleculeAutoCorrect() {
  python::scope module_scope = python::scope();

  python::def(
    "DefaultMoleculePerturber", MoleculeAutoCorrect::DefaultMoleculePerturber);

  boost::python::pointer_wrapper<const MolecularConstraints*>
    null_constraints (nullptr);

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
      &MoleculeAutoCorrect::Settings::max_tree_depth)
    .def_readwrite("n_solutions", 
      &MoleculeAutoCorrect::Settings::n_solutions)
    .def_readwrite("n_top_solutions", 
      &MoleculeAutoCorrect::Settings::n_top_solutions);

  // Define Policy submodule.
  std::string policy_module_name = 
    python::extract<std::string>(python::scope().attr("__name__") + ".Policy");
  python::object policy_module (
    python::handle<>(python::borrowed(
      PyImport_AddModule(policy_module_name.c_str()))));
  python::scope().attr("Policy") = policy_module;
  python::scope policy_scope = policy_module;

  python::enum_<MoleculeAutoCorrect::Policy::Type>("Type")
    .value("BFS", 
      MoleculeAutoCorrect::Policy::Type::BFS)
    .value("Familiarity1", 
      MoleculeAutoCorrect::Policy::Type::Familiarity1)
    .value("Familiarity2", 
      MoleculeAutoCorrect::Policy::Type::Familiarity2)
    .value("DistanceNormalizedFamiliarity", 
      MoleculeAutoCorrect::Policy::Type::DistanceNormalizedFamiliarity)
    .value("Astar", 
      MoleculeAutoCorrect::Policy::Type::Astar)
    .value("UCT", 
      MoleculeAutoCorrect::Policy::Type::UCT)
    .value("MLR", 
      MoleculeAutoCorrect::Policy::Type::MLR);

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

  python::def<
    MoleculeAutoCorrect::Result(
      const RDKit::ROMol&,
      const MoleculeAutoCorrect::Settings&,
      MoleculeAutoCorrect::Policy::Type)>
    ("AutoCorrectMolecule", AutoCorrectMolecule, (
    python::arg("molecule"),
    python::arg("settings"),
    python::arg("selection_policy_type")));

  python::def<
    MoleculeAutoCorrect::Result(
      const RDKit::ROMol&,
      const MoleculeAutoCorrect::Settings&,
      const python::object& objective)>
    ("AutoCorrectMolecule", AutoCorrectMolecule, (
    python::arg("molecule"),
    python::arg("settings"),
    python::arg("objective")));
};

#endif //_PY_MOLECULE_AUTO_CORRECT_HPP_
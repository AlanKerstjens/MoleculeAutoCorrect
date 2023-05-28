#include "MoleculeAutoCorrect.hpp"
#include "MolecularPerturbationUtils.hpp"
#include <chrono>
#include <iostream>

int main(int argc, char* argv[]) {
  if (argc < 3) {
		return 1;
	};

  std::string dictionary_path (argv[1]);
	std::string input_smiles (argv[2]);
  std::size_t n_solutions = 1;
  if (argc > 3) {
    n_solutions = std::stoi(argv[3]);
  };
  bool only_decorate = false;
  if (argc > 4) {
    only_decorate = std::stoi(argv[4]);
  };
  std::size_t n_top_solutions = n_solutions > 10 ? n_solutions : 10;

  ChemicalDictionary dictionary (dictionary_path);

  MoleculePerturber perturber = MoleculeAutoCorrect::DefaultMoleculePerturber();
  if (only_decorate) {
    perturber.SetDecorationSettings();
  };

  MoleculeAutoCorrect::Settings settings (&perturber, &dictionary);
  settings.tree_policy_type = 
    MoleculeAutoCorrect::Settings::TreePolicyType::MLR;

  RDKit::RWMOL_SPTR molecule;
  try {
    molecule = UnsanitizedMoleculeFromSMILES(input_smiles);
    PartialSanitization(*molecule);
  } catch (const RDKit::SmilesParseException& exception) {
    std::cout << exception.what() << std::endl;
    return 1;
  };

  auto begin_time = std::chrono::high_resolution_clock::now();
  MoleculeAutoCorrect::Result result = AutoCorrectMolecule(
    *molecule, settings, n_solutions, n_top_solutions);
  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>
    (end_time - begin_time);

  std::cout << "Search terminated after " << result.n_expansions
            << " expansions (" << elapsed_time.count() << "ms)\n";
  std::cout << "Displaying top " << result.top_discoveries.size()
            << " discoveries" << std::endl;
  for (const MoleculeAutoCorrect::MoleculeDiscovery& discovery :
    result.top_discoveries) {
    std::cout << "Molecule: " << UnsanitizedMoleculeToSMILES(discovery.molecule)
              << ", Familiarity: " << discovery.familiarity
              << ", Discovery order: " << discovery.order
              << ", Discovery depth: " << discovery.depth << std::endl;
  };

  return 0;
};

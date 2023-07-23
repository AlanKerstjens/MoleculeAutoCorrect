#include "MoleculeAutoCorrect.hpp"
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <fstream>

std::pair<std::vector<AtomKey::Error>, std::vector<BondKey::Error>>
MolecularKeyErrors(
	const std::vector<AtomKey>& atom_keys,
	const std::vector<BondKey>& bond_keys,
	const ChemicalDictionary& dictionary) {
	std::size_t n_atoms = atom_keys.size(), n_bonds = bond_keys.size();
	std::vector<AtomKey::Error> atom_key_errors (n_atoms);
	std::vector<BondKey::Error> bond_key_errors (n_bonds);
	for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
		atom_key_errors[atom_idx] = dictionary.AtomKeyError(atom_keys[atom_idx]);
	};
	for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
		bond_key_errors[bond_idx] = dictionary.BondKeyError(bond_keys[bond_idx]);
	};
	return {std::move(atom_key_errors), std::move(bond_key_errors)};
};

std::pair<std::vector<unsigned>, std::vector<unsigned>>
MolecularErroneousness(
	const std::vector<AtomKey::Error>& atom_key_errors,
	const std::vector<BondKey::Error>& bond_key_errors,
	const std::vector<CircularAtomicEnvironment>& foreign_environments){
	std::size_t n_atoms = atom_key_errors.size(), n_bonds = bond_key_errors.size();
	std::vector<unsigned> atom_erroneousness (n_atoms);
	std::vector<unsigned> bond_erroneousness (n_bonds);
	for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
		atom_erroneousness[atom_idx] = static_cast<
			std::underlying_type_t<AtomKey::Error>>(atom_key_errors[atom_idx]);
	};
	for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
		bond_erroneousness[bond_idx] = static_cast<
			std::underlying_type_t<BondKey::Error>>(bond_key_errors[bond_idx]);
	};
	auto [atom_overlaps, bond_overlaps] =
		CircularEnvironmentOverlap(foreign_environments);
	for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
		atom_erroneousness[atom_idx] += atom_overlaps[atom_idx];
	};
	for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
		bond_erroneousness[bond_idx] += bond_overlaps[bond_idx];
	};
	return {std::move(atom_erroneousness), std::move(bond_erroneousness)};
};

int main(int argc, char* argv[]) {
	if (argc != 4) {
		return 1;
	};

  std::string dictionary_path (argv[1]);
	std::string input_smiles (argv[2]);
	std::string output_svg_path (argv[3]);
	bool verbose = true;

  ChemicalDictionary dictionary (dictionary_path);

	RDKit::RWMOL_SPTR molecule;
  try {
    molecule = UnsanitizedMoleculeFromSMILES(input_smiles);
  } catch (const RDKit::SmilesParseException& exception) {
    std::cout << exception.what() << std::endl;
    return 1;
  };

	std::cout << "Input molecule: " << RDKit::MolToSmiles(*molecule) << std::endl;

	std::size_t n_atoms = molecule->getNumAtoms();
	std::size_t n_bonds = molecule->getNumBonds();

	MolecularKeys molecular_keys (*molecule);

	auto [atom_key_errors, bond_key_errors] = MolecularKeyErrors(
		molecular_keys.atom_keys, molecular_keys.bond_keys, dictionary);

	unsigned n_atom_key_errors = std::count_if(
		atom_key_errors.cbegin(), atom_key_errors.cend(),
		[] (AtomKey::Error error) {
			return error != AtomKey::Error::NONE;
		});
	if (n_atom_key_errors > 0) {
		std::cout << "Molecule has " << n_atom_key_errors
							<< " foreign atom keys" << std::endl;
		if (verbose) {
			for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
				AtomKey::Error error = atom_key_errors[atom_idx];
				if (error == AtomKey::Error::NONE) {
					continue;
				};
				const AtomKey& atom_key = molecular_keys.atom_keys[atom_idx];
				std::cout << "Atom " << atom_idx << ": ";
				switch (error) {
					case AtomKey::Error::D:
						std::cout << "D " << atom_key.d();
						break;
					case AtomKey::Error::DV:
						std::cout << "DV " << atom_key.dv();
						break;
					case AtomKey::Error::DVZ:
						std::cout << "DVZ " << atom_key.dvz();
						break;
					case AtomKey::Error::DVZQ:
						std::cout << "DVZQ " << atom_key.dvzq();
						break;
					case AtomKey::Error::DVZQH:
						std::cout << "DVZQH " << atom_key.str();
				};
				std::cout << std::endl;
			};
		};
	};

	unsigned n_bond_key_errors = std::count_if(
		bond_key_errors.cbegin(), bond_key_errors.cend(),
		[] (BondKey::Error error) {
			return error != BondKey::Error::NONE;
		});
	if (n_bond_key_errors > 0) {
		std::cout << "Molecule has " << n_bond_key_errors
							<< " foreign bond keys" << std::endl;
		if (verbose) {
			for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
				BondKey::Error error = bond_key_errors[bond_idx];
				if (error == BondKey::Error::NONE) {
					continue;
				};
				const BondKey& bond_key = molecular_keys.bond_keys[bond_idx];
				std::cout << "Bond " << bond_idx << ": ";
				switch (error) {
					case BondKey::Error::K1K2:
						std::cout << bond_key.k1k2();
						break;
					case BondKey::Error::K1K2B:
						std::cout << bond_key.str();
				};
				std::cout << std::endl;
			};
		};
	};

	std::vector<CircularAtomicEnvironment> foreign_environments =
    MoleculeAutoCorrect::ForeignEnvironments(*molecule, dictionary);

  if (foreign_environments.empty()) {
    std::cout << "Molecule has no foreign circular atomic environments."
							<< std::endl;
    return 0;
  } else {
    std::cout << "Molecule has " << foreign_environments.size()
              << " foreign circular atomic environments." << std::endl;
		if (verbose) {
			for (const CircularAtomicEnvironment& environment : foreign_environments) {
				std::cout << "Environment of atom " << environment.root_atom->getIdx()
									<< ": " << environment.SMILES() << std::endl;
			};
		};
  };

  auto [atom_errors, bond_errors] = MolecularErroneousness(
		atom_key_errors, bond_key_errors, foreign_environments);
	unsigned min_error = std::numeric_limits<unsigned>::max();
	unsigned max_error = 0;
  for (unsigned error : atom_errors) {
    if (error < min_error) {
      min_error = error;
    };
    if (error > max_error) {
      max_error = error;
    };
  };
  for (unsigned error : bond_errors) {
    if (error < min_error) {
      min_error = error;
    };
    if (error > max_error) {
      max_error = error;
    };
  };

  RDKit::DrawColour base_colour (1.0, 0.5, 0.5, 1.0); // RGBA
  double alpha_slope = (1.0 - 0.1) / (max_error - min_error);

  std::vector<int> highlighted_atom_idxs;
  std::map<int, RDKit::DrawColour> atom_highlight_colours;
  for (std::size_t atom_idx = 0; atom_idx < n_atoms; ++atom_idx) {
    unsigned error = atom_errors[atom_idx];
    if (!error) {
      continue;
    };
    highlighted_atom_idxs.push_back(atom_idx);
    RDKit::DrawColour colour (base_colour);
    colour.a = 1.0 - (max_error - error) * alpha_slope;
    atom_highlight_colours.emplace(atom_idx, std::move(colour));
  };

  std::vector<int> highlighted_bond_idxs;
  std::map<int, RDKit::DrawColour> bond_highlight_colours;
  for (std::size_t bond_idx = 0; bond_idx < n_bonds; ++bond_idx) {
    unsigned error = bond_errors[bond_idx];
    if (!error) {
      continue;
    };
    highlighted_bond_idxs.push_back(bond_idx);
    RDKit::DrawColour colour (base_colour);
    colour.a = 1.0 - (max_error - error) * alpha_slope;
    bond_highlight_colours.emplace(bond_idx, std::move(colour));
  };

  std::ofstream svg_stream (output_svg_path);
  RDDepict::compute2DCoords(*molecule);
  RDKit::MolDraw2DSVG drawer(600, 600, svg_stream);
	if (verbose) {
		RDKit::MolDrawOptions& draw_options = drawer.drawOptions();
		draw_options.addAtomIndices = true;
		draw_options.addBondIndices = true;
	};
  drawer.drawMolecule(
    *molecule,
    "",
    &highlighted_atom_idxs,
    &highlighted_bond_idxs,
    &atom_highlight_colours,
    &bond_highlight_colours);
  drawer.finishDrawing();
  svg_stream.close();

  std::cout << "Errors were highlighted in " << output_svg_path << std::endl;

	return 0;
};

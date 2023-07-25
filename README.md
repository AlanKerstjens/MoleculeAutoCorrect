# About

Spell checker for your molecular graphs. A virtual library of reference correct molecules is used to build a dictionary of allowed chemical features. The chemical features of input molecules are compared against this dictionary. If any invalid features are present the molecule is modified in a controlled way to find a closely related valid molecule.

# Installation

## Installation from source

### Prerequisites

Ensure the following dependencies are installed:

* [RDKit](https://rdkit.org/)
* [Molpert](https://github.com/AlanKerstjens/Molpert)
* [Boost](https://www.boost.org/). You already have this if you installed the RDKit. If you'd like to build the Python bindings make sure Boost.Python is installed.
* [CMake](https://cmake.org/)

### Instructions

The following instructions are for GNU+Linux. For alternative operating systems you'll have to adapt these commands slightly.

```shell
git clone https://github.com/AlanKerstjens/MoleculeAutoCorrect.git
export MOLECULE_AUTO_CORRECT="$(pwd)/MoleculeAutoCorrect"
mkdir ${MOLECULE_AUTO_CORRECT}/build && cd ${MOLECULE_AUTO_CORRECT}/build
```

You need to point CMake to your Molpert installation. Assuming it's installed at `${MOLPERT}`:

```shell
cmake -DMolpert_INCLUDE_DIRS=${MOLPERT}/source ..
make install
```

To be able to import the library from Python add `${MOLECULE_AUTO_CORRECT}/lib` to your `${PYTHONPATH}`. Consider doing so in your `bash_profile` file. Otherwise you'll have to manually extend `${PYTHONPATH}` everytime you open a new shell.

```shell
export PYTHONPATH="${PYTHONPATH}:${MOLECULE_AUTO_CORRECT}/lib"
```

### Troubleshooting

CMake will try to find the rest of the dependencies for you. To avoid problems ensure you build the software with the same Boost and Python versions that you used to build Molpert and the RDKit. If CMake finds a different Boost or Python installation you'll need to point it to the correct one, as described [here](https://cmake.org/cmake/help/latest/module/FindBoost.html) and [here](https://cmake.org/cmake/help/latest/module/FindPython.html).

CMake will search for the RDKit in the active Anaconda environment (if you have one) and at `${RDBASE}` if set. If neither of these are the case you need to specify the path to the RDKit yourself. Replace the above CMake command with the one below, substituting the `<placeholder/path>` with your paths.

```shell
cmake -DRDKit_ROOT=<path/to/rdkit> -DMolpert_INCLUDE_DIRS=<path/to/molpert> ..
```

# Quick start

Get your hands on a virtual library of molecules you would like to use as reference of correct chemistry (here `chembl.smi`). Then use this library to create a dictionary of chemical features (here `chembl.dict`). You can specify the radius of circular atomic environments as the last argument (here `1`).

```shell
${MOLECULE_AUTO_CORRECT}/bin/MakeChemicalDictionary chembl.smi chembl.dict 1
```

Given the SMILES string of a molecule (here `[OH2+]C1C([O-])[C]11=C(C#P)C2=S(=C1)C=NC=N2`) you can inspect if it has any issues:

```shell
${MOLECULE_AUTO_CORRECT}/bin/HighlightMoleculeErrors chembl.dict "[OH2+]C1C([O-])[C]11=C(C#P)C2=S(=C1)C=NC=N2" molecule_errors.svg
```

If it has issues you can proceed to try correcting them: 

```shell
python ${MOLECULE_AUTO_CORRECT}/AutoCorrectMolecule.py chembl.dict "[OH2+]C1C([O-])[C]11=C(C#P)C2=S(=C1)C=NC=N2"
```

You can experiment with different settings, including tree policies. Access the `--help` for more information.

```shell
python ${MOLECULE_AUTO_CORRECT}/AutoCorrectMolecule.py --help
```
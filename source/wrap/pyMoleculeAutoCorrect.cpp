#include "pyChemicalDictionary.hpp"
#include "pyMoleculeAutoCorrect.hpp"
#include <boost/python.hpp>

namespace python = boost::python;

BOOST_PYTHON_MODULE(MoleculeAutoCorrect) {
  WrapChemicalDictionary();
  WrapMoleculeAutoCorrect();
};

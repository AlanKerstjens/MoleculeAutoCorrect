cmake_minimum_required(VERSION 3.18)

# Define the project.
project(MoleculeAutoCorrect
  LANGUAGES CXX
  DESCRIPTION "Library for molecular graph correction")

# Specify the C++ standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define the directories where the source code can be found and where the
# compiled libraries/executables should be stored.
set(MOLECULEAUTOCORRECT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/source)
set(MOLECULEAUTOCORRECT_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(MOLECULEAUTOCORRECT_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

# Include the library directory in the default RPATH of the targets.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH ${MOLECULEAUTOCORRECT_LIB_DIR})

# Search for the RDKit.
# If found this will also retrieve the Boost version used to build the RDKit.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(RDKit REQUIRED)

# Check that Molpert was found.
# It's the responsibility of the user to specify its location.
if(NOT DEFINED Molpert_INCLUDE_DIRS OR NOT IS_DIRECTORY ${Molpert_INCLUDE_DIRS})
  message(FATAL_ERROR "Variable Molpert_INCLUDE_DIRS (=\"${Molpert_INCLUDE_DIRS}\") is not set or doesn't specify a valid directory.")
else()
  set(Molpert_FOUND True)
endif()

# The RDKit may have been compiled without Boost.Serialization support.
# Ensure we have Boost.Serialization.
find_package(Boost REQUIRED COMPONENTS serialization)

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${RDKit_INCLUDE_DIRS}
  ${Molpert_INCLUDE_DIRS}
  ${MOLECULEAUTOCORRECT_SOURCE_DIR})

link_directories(
  ${Boost_LIBRARY_DIRS}
  ${RDKit_LIBRARY_DIRS})

# Define the library targets.
add_library(MoleculeAutoCorrect INTERFACE)

target_link_libraries(MoleculeAutoCorrect INTERFACE
  Boost::serialization
  RDKitRDGeneral
  RDKitGraphMol
  RDKitFingerprints
  RDKitSmilesParse)

option(BUILD_PYTHON_BINDINGS "Enable Python binding compilation" ON)
if (BUILD_PYTHON_BINDINGS)
  find_package(Python REQUIRED COMPONENTS Interpreter Development)
  find_package(Boost REQUIRED COMPONENTS python)

  add_library(pyMoleculeAutoCorrect MODULE
    ${MOLECULEAUTOCORRECT_SOURCE_DIR}/wrap/pyMoleculeAutoCorrect.cpp)

  target_include_directories(pyMoleculeAutoCorrect PUBLIC
    ${Molpert_INCLUDE_DIRS}/wrap
    ${MOLECULEAUTOCORRECT_SOURCE_DIR}/wrap
    ${Python_INCLUDE_DIRS})

  target_link_directories(pyMoleculeAutoCorrect PUBLIC
    ${Python_LIBRARY_DIRS})

  target_link_libraries(pyMoleculeAutoCorrect PUBLIC
    ${Python_LIBRARIES}
    Boost::python
    MoleculeAutoCorrect)

  set_target_properties(pyMoleculeAutoCorrect PROPERTIES
    PREFIX ""
    IMPORT_PREFIX ""
    OUTPUT_NAME "MoleculeAutoCorrect")

  install(
    TARGETS
      pyMoleculeAutoCorrect
    DESTINATION ${MOLECULEAUTOCORRECT_LIB_DIR})
endif()

# Define the executable targets
add_executable(AutoCorrectMolecule
  ${MOLECULEAUTOCORRECT_SOURCE_DIR}/AutoCorrectMolecule.cpp)

target_link_libraries(AutoCorrectMolecule PRIVATE
  MoleculeAutoCorrect)

add_executable(MakeChemicalDictionary
  ${MOLECULEAUTOCORRECT_SOURCE_DIR}/MakeChemicalDictionary.cpp)

target_link_libraries(MakeChemicalDictionary PRIVATE
  MoleculeAutoCorrect
  RDKitFileParsers)

add_executable(HighlightMoleculeErrors
  ${MOLECULEAUTOCORRECT_SOURCE_DIR}/HighlightMoleculeErrors.cpp)

target_link_libraries(HighlightMoleculeErrors PRIVATE
  MoleculeAutoCorrect
  RDKitDepictor
  RDKitMolDraw2D)

# Specify where the executables should be installed.
install(
  TARGETS
    AutoCorrectMolecule
    MakeChemicalDictionary
    HighlightMoleculeErrors
  DESTINATION ${MOLECULEAUTOCORRECT_BIN_DIR})
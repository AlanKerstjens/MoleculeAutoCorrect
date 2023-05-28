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
set(MOLECULEAUTOCORRECT_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(MOLECULEAUTOCORRECT_SOURCE_DIR ${MOLECULEAUTOCORRECT_ROOT_DIR}/source)
set(MOLECULEAUTOCORRECT_BIN_DIR ${MOLECULEAUTOCORRECT_ROOT_DIR}/bin)
set(MOLECULEAUTOCORRECT_LIB_DIR ${MOLECULEAUTOCORRECT_ROOT_DIR}/lib)

# Include the library directory in the default RPATH of the targets.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH ${MOLECULEAUTOCORRECT_LIB_DIR})

# Check that the RDKit was found.
# It's the responsibility of the user to specify its location.
if(NOT DEFINED RDKit_INCLUDE_DIRS OR NOT IS_DIRECTORY ${RDKit_INCLUDE_DIRS})
  message(FATAL_ERROR "Variable RDKit_INCLUDE_DIRS (=\"${RDKit_INCLUDE_DIRS}\") is not set or doesn't specify a valid directory.")
elseif(NOT DEFINED RDKit_LIBRARY_DIRS OR NOT IS_DIRECTORY ${RDKit_LIBRARY_DIRS})
  message(FATAL_ERROR "Variable RDKit_LIBRARY_DIRS (=\"${RDKit_LIBRARY_DIRS}\") is not set or doesn't specify a valid directory.")
else()
  set(RDKit_FOUND True)
endif()

# Check that Molpert was found.
# It's the responsibility of the user to specify its location.
if(NOT DEFINED Molpert_INCLUDE_DIRS OR NOT IS_DIRECTORY ${Molpert_INCLUDE_DIRS})
  message(FATAL_ERROR "Variable Molpert_INCLUDE_DIRS (=\"${Molpert_INCLUDE_DIRS}\") is not set or doesn't specify a valid directory.")
else()
  set(Molpert_FOUND True)
endif()

# Search for the rest of the dependencies.
find_package(Boost REQUIRED COMPONENTS serialization)

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${RDKit_INCLUDE_DIRS}
  ${Molpert_INCLUDE_DIRS}
  ${Molpert_INCLUDE_DIRS}/wrap
  ${MOLECULEAUTOCORRECT_SOURCE_DIR}
  ${MOLECULEAUTOCORRECT_SOURCE_DIR}/wrap)

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

# Specify where the executables should be installed.
install(
  TARGETS
    AutoCorrectMolecule
    MakeChemicalDictionary
  DESTINATION ${MOLECULEAUTOCORRECT_BIN_DIR})
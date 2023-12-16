import os
import time
import argparse
from rdkit import Chem
import molpert as mpt
import MoleculeAutoCorrect as mac

def ParseArgs():
    parser = argparse.ArgumentParser(
        description="Corrects an input molecular graph",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dictionary", type=str,
        help="Path to chemical dictionary.")
    parser.add_argument("smiles", type=str,
        help="Input molecule SMILES string.")
    parser.add_argument("-p", "--policy", type=str, default="MLR",
        choices=["BFS", "Familiarity1", "Familiarity2",
                 "DistanceNormalizedFamiliarity", "UCT", "Astar", "MLR"],
        help="Tree search vertex selection policy type.")
    parser.add_argument("-s", "--max_tree_size", type=int, default=1000,
        help="Maximum tree size.")
    parser.add_argument("-d", "--max_tree_depth", type=int, default=10,
        help="Maximum tree depth.")
    parser.add_argument("-n", "--n_solutions", type=int, default=1,
        help="Number of requested solutions.")
    parser.add_argument("-do", "--decorate_only", action="store_true",
        help="Flag to only enable decorative perturbations.")
    parser.add_argument("-da", "--disable_atom_insertions", action="store_true",
        help="Flag to disable atom insertions.")
    args = parser.parse_args()
    if not os.path.isfile(args.dictionary):
        raise parser.error("dictionary isn't a file.")
    return args

def Main():
    args = ParseArgs()

    dictionary = mac.ChemicalDictionary(args.dictionary)

    perturber = mac.DefaultMoleculePerturber()
    if args.decorate_only:
        perturber.SetDecorationSettings()

    settings = mac.Settings(perturber, dictionary)
    settings.max_tree_size = args.max_tree_size
    settings.max_tree_depth = args.max_tree_depth
    settings.n_solutions = args.n_solutions
    settings.n_top_solutions = args.n_solutions
    if args.disable_atom_insertions:
        settings.attempt_degree_correction_with_insertions = False
        settings.attempt_bond_correction_with_atom_insertions = False
        settings.attempt_environment_correction_with_atom_insertions = False

    policy_types = {
        "BFS": mac.Policy.Type.BFS,
        "Familiarity1": mac.Policy.Type.Familiarity1,
        "Familiarity2": mac.Policy.Type.Familiarity2,
        "DistanceNormalizedFamiliarity": mac.Policy.Type.DistanceNormalizedFamiliarity,
        "UCT": mac.Policy.Type.UCT,
        "Astar": mac.Policy.Type.Astar,
        "MLR": mac.Policy.Type.MLR}

    molecule = Chem.MolFromSmiles(args.smiles, sanitize=False)

    start_time = time.perf_counter_ns()
    mpt.PartialSanitization(molecule)
    result = mac.AutoCorrectMolecule(
        molecule=molecule,
        settings=settings,
        selection_policy_type=policy_types[args.policy])
    end_time = time.perf_counter_ns()
    
    print(f"Search terminated after {result.n_expansions} expansions ({(end_time - start_time)/1e6:.0f} ms)")
    for discovery in result.top_discoveries:
        print((f"Molecule: {Chem.MolToSmiles(discovery.molecule)}, "
               f"Familiarity: {discovery.familiarity}, "
               f"Order: {discovery.order}, "
               f"Depth: {discovery.depth}"))

if __name__ == "__main__":
    Main()
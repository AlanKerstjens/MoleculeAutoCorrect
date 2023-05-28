import os
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
    parser.add_argument("-p", "--policy", type=str, default="Familiarity",
        choices=["Familiarity", "BFS", "DistanceNormalizedFamiliarity",
                 "UCT", "Astar", "MLR"],
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
        help="Flag to disable atom insertions in bond and environment correction.")
    args = parser.parse_args()
    if not os.path.isfile(args.dictionary):
        raise parser.error("dictionary isn't a file.")
    return args

def TreePolicyType(policy_name):
    if policy_name == "Familiarity":
        return mac.Settings.TreePolicyType.Familiarity
    elif policy_name == "BFS":
        return mac.Settings.TreePolicyType.BFS
    elif policy_name == "DistanceNormalizedFamiliarity":
        return mac.Settings.TreePolicyType.DistanceNormalizedFamiliarity
    elif policy_name == "UCT":
        return mac.Settings.TreePolicyType.UCT
    elif policy_name == "Astar":
        return mac.Settings.TreePolicyType.Astar
    elif policy_name == "MLR":
        return mac.Settings.TreePolicyType.MLR
    raise ValueError("Unknown policy name")

def Main():
    args = ParseArgs()

    dictionary = mac.ChemicalDictionary(args.dictionary)

    perturber = mac.DefaultMoleculePerturber()
    if args.decorate_only:
        perturber.SetDecorationSettings()

    settings = mac.Settings(perturber, dictionary)
    settings.tree_policy_type = TreePolicyType(args.policy)
    settings.max_tree_size = args.max_tree_size
    settings.max_tree_depth = args.max_tree_depth
    if args.disable_atom_insertions:
        settings.attempt_bond_correction_with_atom_insertions = False
        settings.attempt_environment_correction_with_atom_insertions = False

    molecule = Chem.MolFromSmiles(args.smiles, sanitize=False)
    mpt.PartialSanitization(molecule)

    result = mac.AutoCorrectMolecule(molecule, settings, 
        n_solutions=args.n_solutions, n_top_solutions=args.n_solutions)
    for discovery in result.top_discoveries:
        print((f"Molecule: {Chem.MolToSmiles(discovery.molecule)}, "
               f"Familiarity: {discovery.familiarity}, "
               f"Order: {discovery.order}, "
               f"Depth: {discovery.depth}"))

if __name__ == "__main__":
    Main()
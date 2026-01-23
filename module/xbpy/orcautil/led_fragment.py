from scipy.spatial import cKDTree
from rdkit.Chem import AllChem as Chem
from ..rdutil import position
from ..rdutil import get_connected_component_indices

def fragment_molecule(mol, center_atom_index):
    """
    Assumes the molecule is already proximity bonded.
    """
    all_components = get_connected_component_indices(mol)
    center_component_indices = None
    for component_indices in all_components:
        if center_atom_index in component_indices:
            center_component_indices = component_indices
            break
    if center_component_indices is None:
        raise ValueError(f"Center atom index {center_atom_index} not found in molecule")

    surrounding_indices = get_surrounding_indices(mol, center_component_indices)

    # extend surrounding indices to next C-C bond.
    all_used_components = []
    for component_indices in all_components:
        current_selected = set(component_indices).intersection(surrounding_indices)
        current_selected = extend_connected_indices(mol, current_selected, component_indices)
        all_used_components.append(current_selected)

    return all_used_components

def get_surrounding_indices(mol, center_indices, radius = 5):
    positions = position(mol)
    tree = cKDTree(positions)
    surrounding_indices = tree.query_ball_point(positions[center_indices], radius)
    surrounding_indices = [index for sublist in surrounding_indices for index in sublist]
    return surrounding_indices

def extend_connected_indices(mol, indices, component_indices):
    """ Extends the indices within the connected component of the given indices to the next C-C bond."""
    to_check = list(indices)
    seen = set(indices)
    while to_check:
        index = to_check.pop()
        for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            single_bond = mol.GetBondBetweenAtoms(index, neighbor_index).GetBondType() == Chem.rdchem.BondType.SINGLE
            both_carbon = mol.GetAtomWithIdx(index).GetSymbol() == "C" and mol.GetAtomWithIdx(neighbor_index).GetSymbol() == "C"
            seen.add(neighbor_index)
            if not(single_bond and both_carbon):
                to_check.append(neighbor_index)
    return list(seen)

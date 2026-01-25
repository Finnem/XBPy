from scipy.spatial import cKDTree
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdmolops
from ..rdutil import position
from ..rdutil import get_connected_component_indices
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def fragment_molecule(mol, center_atom_index, verbose = False):
    """
    Assumes the molecule is already proximity bonded.
    """
    all_components = get_connected_component_indices(mol)
    center_component_indices = None
    for component_indices in all_components:
        if center_atom_index in component_indices:
            center_component_indices = component_indices
            break
    if verbose:
        logger.info(f"Found center component {center_component_indices}")
    if center_component_indices is None:
        raise ValueError(f"Center atom index {center_atom_index} not found in molecule")

    surrounding_indices = get_surrounding_indices(mol, center_component_indices)

    # extend surrounding indices to next C-C bond.
    all_used_components = []
    all_to_replace = []
    ring_info = None
    try:
        rdmolops.FastFindRings(mol)
        ring_info = mol.GetRingInfo()
    except Exception:
        ring_info = None

    for component_indices in all_components:
        current_selected = set(component_indices).intersection(surrounding_indices)
        current_selected, to_replace = extend_connected_indices(
            mol, current_selected, component_indices, ring_info=ring_info, max_ring_size=12
        )
        if len(current_selected) != 0:
            all_used_components.append(current_selected)
            if verbose:
                logger.info(f"Extended component {len(component_indices)} to {len(current_selected)}")
            all_to_replace.extend(to_replace)
    return all_used_components, list(set(all_to_replace))

def get_surrounding_indices(mol, center_indices, radius = 5):
    positions = position(mol)
    tree = cKDTree(positions)
    surrounding_indices = tree.query_ball_point(positions[center_indices], radius)
    surrounding_indices = [index for sublist in surrounding_indices for index in sublist]
    return surrounding_indices

def extend_connected_indices(mol, indices, component_indices, ring_info=None, max_ring_size=12):
    """ Extends the indices within the connected component of the given indices to the next C-C bond."""
    to_check = list(indices)
    seen = set(indices)
    keep = set(indices)
    component_set = set(component_indices)
    small_ring_atoms = set()
    if ring_info is not None:
        try:
            for idx in component_set:
                for size in range(3, max_ring_size + 1):
                    if ring_info.IsAtomInRingOfSize(int(idx), int(size)):
                        small_ring_atoms.add(idx)
                        break
        except Exception:
            small_ring_atoms = set()
    # heavy atom check
    while to_check:
        index = to_check.pop()
        index_atom = mol.GetAtomWithIdx(index)
        for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor_index not in component_set:
                continue
            neighbor_symbol = neighbor.GetSymbol()
            index_symbol = index_atom.GetSymbol()
            single_bond = mol.GetBondBetweenAtoms(index, neighbor_index).GetBondType() == Chem.rdchem.BondType.SINGLE
            both_carbon = (index_symbol == "C") and (neighbor_symbol == "C")
            # do not cut rings
            if not (neighbor_index in seen) and (not(single_bond and both_carbon and not(neighbor_index in small_ring_atoms) and not(index in small_ring_atoms))):
                to_check.append(neighbor_index)
            seen.add(neighbor_index)
            if both_carbon:
                # check if neighbor is sp3, we proxy by having 4 neighbors since we know its carbon
                if len(neighbor.GetNeighbors()) == 4:
                    keep.add(neighbor_index)
            else:
                keep.add(neighbor_index)
    # add any hydrogen neighbors
    new_seen = set()
    for index in keep:
        for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 1:
                new_seen.add(neighbor_index)
    
    # determine atoms to replace with hydrogens
    to_replace = set()
    for index in keep:
        for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor_index not in keep:
                if neighbor.GetAtomicNum() != 1:
                    to_replace.add((index, neighbor_index))
    keep.update(new_seen)
    return list(keep), list(to_replace)

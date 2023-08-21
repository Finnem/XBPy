import numpy as np
import rdkit.Chem as Chem
from rdkit.Chem import rdmolops
from scipy.stats import rankdata
periodic_table = Chem.GetPeriodicTable()
from collections import defaultdict
import logging

from ..rdutil import position

def unique_index(mol):
    connected_components = _connected_components(mol)


    original_scores = []
    orders = []
    # determine scores for each connected component
    for connected_component in connected_components:
        original_scores.append(morgan_prop(mol, connected_component))
        # resolve ambiguities
        order = rankdata(-np.round(original_scores[-1], 4), method="min") - 1
        orders.append(resolve_ambiguity([mol.GetAtomWithIdx(i) for i in connected_component], order))

    # sort connected components based on scores. For some reason, lexsort sorts from last to first
    max_len = max([len(s) for s in original_scores])
    padded_scores = np.array([np.pad(s, (0, max_len - len(s)), mode="constant", constant_values=0) for s in original_scores])
    component_order = np.lexsort(np.array([list(reversed(s)) for s in padded_scores]).T)

    new_index = np.zeros(mol.GetNumAtoms(), dtype=int)
    # determine final order
    cur_index_offset = 0
    for i in component_order:
        new_index[connected_components[i]] = orders[i] + cur_index_offset
        cur_index_offset += len(connected_components[i])    
    return new_index

def _connected_components(mol):
    considered_indices = set(range(mol.GetNumAtoms()))
    components = []
    while considered_indices:
        component = set()
        queue = [considered_indices.pop()]
        while queue:
            index = queue.pop()
            component.add(index)
            for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
                neighbor_index = neighbor.GetIdx()
                if neighbor_index in considered_indices:
                    queue.append(neighbor_index)
                    considered_indices.remove(neighbor_index)
        components.append(list(component))
    return components

def morgan_prop(mol, indices = None, return_history = False):
    if indices is None:
        indices = list(range(mol.GetNumAtoms()))
    
    atoms = [mol.GetAtomWithIdx(i) for i in indices]

    # initial morgan propagation
    adjacency_matrix = np.array(rdmolops.GetAdjacencyMatrix(mol), dtype=float)[indices, :][:, indices]
    old_err = np.geterr()["invalid"]
    np.seterr(invalid="ignore")
    adjacency_matrix /= adjacency_matrix.sum(axis=1, keepdims=True)
    np.seterr(invalid=old_err)
    adjacency_matrix = np.nan_to_num(adjacency_matrix, nan=0)
    adjacency_matrix += np.eye(adjacency_matrix.shape[0], dtype=int)
    element_vector = np.array([a.GetAtomicNum() for a in atoms], dtype = float)
    last_unique = None
    if return_history: history = [np.zeros(len(atoms)), element_vector]
    while last_unique != (last_unique := len(np.unique(element_vector))):
        element_vector = adjacency_matrix @ element_vector
        if return_history: history.append(element_vector)
    if return_history:
        return np.array(history)
    return element_vector

def get_equivalent_atoms(mol, indices):
    element_vector = morgan_prop(mol, indices)
    order = rankdata(-np.round(element_vector, 4), method="min") - 1
    result = []
    unique_values, unique_counts = np.unique(order, return_counts=True)
    for val in unique_values[unique_counts > 1]:
        result.append(np.where(order == val)[0])
    return result

def resolve_ambiguity(atoms, order, ambiguous_indices = None):

    # resolve ambiguities based on spatial distance to heaviest atoms
    sorted_indices = np.argsort(order)

    def _resolve_ambiguity(ambigous_index, relative_reference_index):
        # ambigous indices are indices of atoms which were assigned the same index up to now
        # relative reference index is the index of the atom which is used as reference for the relative position
        # atoms are then sorted based on their proximity to the reference atom
        
        # if no more atoms are left, go to worst case
        if relative_reference_index < 0:
            return _resolve_ambiguity_worst_case(ambigous_index)
        
        # otherwise determine distances to reference atom
        positions = np.array([position(atoms[int(i)]) for i in ambigous_index])
        target_position = position(atoms[int(sorted_indices[relative_reference_index])])
        distances = np.linalg.norm(np.array(target_position - positions), axis=1)
        results = rankdata(-distances, method="min") - 1
        unique_vals, unique_counts = np.unique(results, return_counts=True)

        # if any ambiguity remains, go to next reference atom
        for val in unique_vals[unique_counts > 1]:
            next_indices = np.where(results == val)[0]
            results[next_indices] += _resolve_ambiguity(ambigous_index[next_indices], relative_reference_index - 1)
        return results

    def _resolve_ambiguity_worst_case(ambigous_index):
        positions = np.array([position(atoms[int(i)]) for i in ambigous_index])
        distances = np.linalg.norm(np.array(positions), axis=1)
        results = np.argsort(-distances) - 1
        return results
    
    if ambiguous_indices is None:
        unique_vals, unique_counts = np.unique(order, return_counts=True)
        for val in unique_vals[unique_counts > 1]:
            next_indices = np.where(order == val)[0]
            order[next_indices] += _resolve_ambiguity(next_indices, len(sorted_indices) - 1)
    
    else:
        return _resolve_ambiguity(ambiguous_indices, len(sorted_indices) - 1)
    return order


def substructure_match(mol1, mol2):
    mol1_history = morgan_prop(mol1, return_history=True)
    mol2_history = morgan_prop(mol2, return_history=True)
    index_map = {}
    indexed_mol2 = set()
    for i in range(min(mol1_history.shape[0], mol2_history.shape[0])):
        # constructing history hash maps:
        history_map1 = defaultdict(lambda: [])
        for index, history_vector in enumerate(mol1_history.T):
            if index in index_map: continue
            history_map1[tuple(history_vector[:-(i + 1)])].append(index)
        history_map2 = defaultdict(lambda: [])
        for index, history_vector in enumerate(mol2_history.T):
            if index in indexed_mol2: continue
            history_map2[tuple(history_vector[:-(i + 1)])].append(index)       
        
        for history_vector, indices in history_map1.items():
            if history_vector in history_map2:
                target_indices = history_map2[history_vector]
                min_length = min(len(target_indices), len(indices))
                if len(target_indices) != len(indices):
                    logging.warning(f"Found ambiguous equivalent atoms matching {indices} in molecule 1 to atoms {target_indices} in molecule 2. Will match {indices[:min_length]} to {target_indices[:min_length]}.")
                used_indices = np.array(indices)[:min_length]
                used_target_indices = np.array(target_indices)[:min_length]
                if len(used_indices) > 1:
                    used_indices = used_indices[resolve_ambiguity(mol1.GetAtoms(), mol1_history[1], used_indices)]
                    used_target_indices = used_target_indices[resolve_ambiguity(mol2.GetAtoms(), mol2_history[2], used_target_indices)]
                for used_index, used_target_index in zip(used_indices, used_target_indices):
                    index_map[used_index] = used_target_index
                    indexed_mol2.add(used_target_index)
    result = np.full(len(mol2.GetAtoms()), -1)
    for index, target in index_map.items():
        result[index] = target
    return result

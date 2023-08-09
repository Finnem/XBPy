import numpy as np
import rdkit.Chem as Chem
from rdkit.Chem import rdmolops
from scipy.stats import rankdata
periodic_table = Chem.GetPeriodicTable()

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

def morgan_prop(mol, indices = None):
    if indices is None:
        indices = list(range(mol.GetNumAtoms()))
    
    atoms = [mol.GetAtomWithIdx(i) for i in indices]

    # initial morgan propagation
    adjacency_matrix = np.array(rdmolops.GetAdjacencyMatrix(mol), dtype=float)[indices, :][:, indices]
    adjacency_matrix /= adjacency_matrix.sum(axis=1, keepdims=True)
    adjacency_matrix += np.eye(adjacency_matrix.shape[0], dtype=int)
    element_vector = np.array([a.GetAtomicNum() for a in atoms], dtype = float)
    last_unique = None
    while last_unique != (last_unique := len(np.unique(element_vector))):
        element_vector = adjacency_matrix @ element_vector
    return element_vector

def get_equivalent_atoms(mol, indices):
    element_vector = morgan_prop(mol, indices)
    order = rankdata(-np.round(element_vector, 4), method="min") - 1
    result = []
    unique_values, unique_counts = np.unique(order, return_counts=True)
    for val in unique_values[unique_counts > 1]:
        result.append(np.where(order == val)[0])
    return result

def resolve_ambiguity(atoms, order):

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
        results = np.argsort(-distances, method="min") - 1
        return results
    
    unique_vals, unique_counts = np.unique(order, return_counts=True)
    for val in unique_vals[unique_counts > 1]:
        next_indices = np.where(order == val)[0]
        order[next_indices] += _resolve_ambiguity(next_indices, len(sorted_indices) - 1)

    return order
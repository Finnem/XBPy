import numpy as np
from ..rdutil import position
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
import json
import os
from rdkit import Chem

def distribute_to_batches(molecules, batch_size, out_folder, format, similarity_threshold=None, suffix = None, start = 0):
    """Distribute molecules to batches, detect similar molecules to create a dependency file and write them to disk.

    Args:
        molecules (list(RDKit.Mol)): List of RDKit molecules.
        batch_size (int): Size of each batch.
        out_path (str): Path to write the batches to.
        format (str): Format to write the batches in. Can be "xyz", "mol" or "pdb".
        similarity_threshold (float): Defaults to None. Threshold for similarity to create a dependency file.
        suffix (str): Defaults to None. Suffix to append to the batch files.

    Returns:
        list(str): List of paths to the written batches.

    """
    raise DeprecationWarning("This function is deprecated. Please use the function in the rdutil.io.write_molecules function.")
    if suffix is None:
        suffix = ""
    if similarity_threshold is None:
        batches = [range(i, min(i+batch_size, len(molecules))) for i in range(0, len(molecules), batch_size)]
    else:
        positions = np.array([position(mol) for mol in molecules])
        mean_positions = np.mean(positions, axis=1)
        tree = cKDTree(mean_positions)
        incidence_matrix = tree.sparse_distance_matrix(tree, max_distance=similarity_threshold)
        cols, rows = incidence_matrix.nonzero()
        pairwise_rmsds = np.array([np.sqrt(np.mean((positions[cols[i]] - positions[rows[i]])**2)) for i in range(len(cols))])
        incidence_matrix = csr_matrix((pairwise_rmsds[pairwise_rmsds < similarity_threshold], np.array([cols, rows])[:,pairwise_rmsds < similarity_threshold]), shape=incidence_matrix.shape)
        n_components, labels = connected_components(incidence_matrix, directed=False)
        # greedily fill up batches
        batches = []
        cur_batch_size = 0
        cur_batch = []
        label_counts = np.bincount(labels)
        for i in range(n_components):
            labelled_indices = np.where(labels == i)[0]
            for j in np.arange(0, label_counts[i], batch_size):
                end = min(j+batch_size, label_counts[i])
                difference = end - j
                cur_batch_size += difference
                if cur_batch_size > batch_size:
                    batches.append(np.concatenate(cur_batch))
                    cur_batch = [labelled_indices[j:end]]
                    cur_batch_size = difference
                else:
                    cur_batch.append(labelled_indices[j:end])
                    cur_batch_size += difference
        batches.append(np.concatenate(cur_batch))

    # write batches to files and create computational dependencies
    paths_to_mols = {}
    os.makedirs(out_folder, exist_ok=True)
    mol_idx = 0
    batch_paths = []
    for i, batch in zip(range(len(batches)), batches):
        out_index = i + start
        batch_path = os.path.join(out_folder, f"batch_{out_index}")
        batch_paths.append(batch_path)
        os.makedirs(batch_path, exist_ok=True)
        for mol_index in batch:
            out_path = os.path.join(batch_path, str(mol_index) + suffix + ".xyz")
            paths_to_mols[mol_index] = out_path
            if format == "mol":
                Chem.MolToMolFile(molecules[mol_index], out_path)
            elif format == "pdb":
                Chem.MolToPDBFile(molecules[mol_index], out_path)
            elif format == "xyz":
                Chem.MolToXYZFile(molecules[mol_index], out_path)
            else:
                raise ValueError("Format not recognized.")
            mol_idx += 1

    # write dependencies for each batch
    if similarity_threshold is not None:
        dense_incidence = np.array(incidence_matrix.todense())
        for i, batch in enumerate(batches):
            with open(os.path.join(out_folder, f"batch_{i}", "dependencies.json"), "w") as f:
                print(i, {str(mol_idx): [(paths_to_mols[j], dense_incidence[mol_idx, j]) for j in dense_incidence[mol_idx].nonzero()[0]] for mol_idx in batch}, f)
                json.dump({str(mol_idx): [(paths_to_mols[j], dense_incidence[mol_idx, j]) for j in dense_incidence[mol_idx].nonzero()[0]] for mol_idx in batch}, f)

    return batch_paths
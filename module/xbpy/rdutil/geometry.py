from rdkit.Chem import rdMolTransforms
from rdkit import Chem
import numpy as np
import logging
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix, coo_matrix

class AtomKDTree():
    def __init__(self, atoms):
        vdw_radii = np.array([Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in atoms])
        elements = np.array([atom.GetAtomicNum() for atom in atoms])
        self.trees = {}
        for vdw_radius, vdw_index in zip(*np.unique(vdw_radii, return_index = True)):
            tree = cKDTree(position(atoms)[vdw_radii == vdw_radius])
            self.trees[elements[vdw_index]] = (vdw_radius, tree)

    def query_ball_point(self, points, radius, element = None):
        indices = [set() for _ in range(len(points))]
        if element == None:
            considered_trees = self.trees.keys()
        else:
            considered_trees = [element]
        for key in considered_trees:
            vdw_radius, tree in self.trees[key]
            query_result = tree.query_ball_point(points, radius + vdw_radius)
            for i, result in enumerate(query_result):
                indices[i].update(result)
        return indices

    def query_ball_tree(self, other, radius, element = None):
        indices = [set() for _ in range(len(self.trees))]
        if element == None:
            considered_trees = self.trees.keys()
        else:
            considered_trees = [element]
        for key in considered_trees:
            if key not in other.trees:
                continue
            vdw_radius, tree = self.trees[key]
            query_result = tree.query_ball_tree(other.trees[key][1], radius + vdw_radius)
            for i, result in enumerate(query_result):
                indices[i].update(result)
        return indices


    def sparse_distance_matrix(self, other, max_distance, element = None):
        if element == None:
            considered_trees = self.trees.keys()
        else:
            considered_trees = [element]
        data = []
        row = []
        col = []
        for key in considered_trees:
            if key not in other.trees:
                continue
            vdw_radius, tree = self.trees[key]
            query_result = tree.sparse_distance_matrix(other.trees[key][1], max_distance + vdw_radius, output_type = 'coo_matrix')
            data.append(query_result.data)
            row.append(query_result.row)
            col.append(query_result.col)
        if sum([len(a) for a in data]) == 0:
            return csr_matrix((0,0))
        return coo_matrix((np.concatenate(data), (np.concatenate(row), np.concatenate(col))))

def detect_clash(atoms1, atoms2):
    """
    Detect if the given atoms are in a clash. A clash is detected if the distance between any two atoms is smaller than the sum of their van-der-waals radii.

    Args:
        atoms1 (list): List of atoms to check for clashes.
        atoms2 (list): List of atoms to check for clashes.

    Returns:
        bool: True if a clash is detected, False otherwise.

    """
    atoms1 = position(atoms1)
    atoms2 = position(atoms2)
    periodic_table = Chem.GetPeriodicTable()
    vdw_radii1 = np.array([periodic_table.GetRvdw(atom.GetAtomicNum()) for atom in atoms1])
    vdw_radii2 = np.array([periodic_table.GetRvdw(atom.GetAtomicNum()) for atom in atoms2])
    distances = np.linalg.norm(atoms1[:, None, :] - atoms2[None, :, :], axis = -1)
    return (distances < vdw_radii1[:, None] + vdw_radii2[None, :]).any()

def positions(*args, **kwargs):
    return position(*args, **kwargs)

def position(atom, conformer = 0, force_seperate_check = False, retry_with_mol_split = True):
    """Return the position of the given atom.

    Args:
        atom (RDKit.Atom): Atom to get the position of. If list will return a list of positions.
        conformer (int): Defaults to 0. Conformer to get the position from.

    Returns:
        np.ndarray: 3xN Position of the given N atoms.

    """
    import numpy as np
    try:
        atom = list(atom)
    except TypeError:
        pass
    if isinstance(atom, list) or isinstance(atom, tuple) or isinstance(atom, np.ndarray):
        try:
            if not (type(conformer) == int):
                conformer = conformer
            else:
                conformer = atom[0].GetOwningMol().GetConformer(conformer)

            if (len(atom) > 10) and not force_seperate_check: #TODO: 10 is currently a magic number, should be replaced by a more sensible value
                # assume its faster to get all positions at once and then slice
                
                molecular_positions = np.array(conformer.GetPositions())
                indices = [a.GetIdx() for a in atom]
                return molecular_positions[indices]
            else:
                return np.array([conformer.GetAtomPosition(a.GetIdx()) for a in atom])
        except (RuntimeError, IndexError):
            # try to split atoms by molecules, compute positions and concatenate
            if retry_with_mol_split:
                atom_by_mol = {}
                for a in atom:
                    mol = a.GetOwningMol()
                    if mol not in atom_by_mol:
                        atom_by_mol[mol] = []
                    atom_by_mol[mol].append(a)
                positions = []
                if len(atom_by_mol) == 1:
                    raise
                for mol, atoms in atom_by_mol.items():
                    positions.extend(position(atoms, retry_with_mol_split = False))
                return np.array(positions)
            else:
                raise
    if isinstance(atom, Chem.rdchem.Mol):
        conformer = atom.GetConformer(conformer)
        return np.array(conformer.GetPositions())
    
    else:
        if not (type(conformer) == int):
            conformer = conformer
        else:
            conformer = atom.GetOwningMol().GetConformer(conformer)
        return np.array(conformer.GetAtomPosition(atom.GetIdx()))
        

def transform(mol, matrix, translation = None):
    """
    Transform the given molecule INPLACE by the given matrix and translation vector.
    
    Args:
        mol (RDKit.Mol): Molecule to transform.
        matrix (np.ndarray): Rotation matrix (3x3) or transformation matrix (4x4).
        translation (np.ndarray): Translation vector (3x1).



    """
    if matrix.shape == (3,3):
        transformation = np.eye(4)
        transformation[:3, :3] = matrix
        if translation is not None:
            transformation[:3, 3] = translation
    if matrix.shape == (4,4):
        transformation = matrix
        if translation is not None:
            logging.warning("Translation vector was explicitly passed to transform, however the 4x4 transformation matrix already implies a translation. The passed translation vector will be ignored.")
    rdMolTransforms.TransformConformer(mol.GetConformer(), transformation)


def check_occlusion(from_atoms, to_atoms, potential_occluders, return_occluders = False, ignore_outside = True, default_radius = 1, occluders_vdw_factor = 1):
    """
    Check if the bond between from_atoms and to_atoms is occluded by any of the atoms in potential_occluders.
    If return_occluders is True, return the occluding atoms. An atom is considered occluding the bond between two other atoms if its van-der-waals radius is larger than its distance to the bond.
    The bonds are assumed to be between the atoms at the same index in from_atoms and to_atoms.

    Args:
        from_atoms (list): List of atoms at the start of the bonds. Alteratively, a position array can be passed. Then the default radius will be used for the atoms.
        to_atoms (list): List of atoms at the end of the bonds. Alteratively, a position array can be passed. Then the default radius will be used for the atoms.
        potential_occluders (list): List of atoms that could occlude the bonds. Alteratively, a position array can be passed. Then the default radius will be used for the atoms.
        return_occluders (bool): Defaults to False. If True, return a boolean array marking each occluding atom for each bond.
        ignore_outside (bool): Defaults to True. If True, only consider atoms that are between the two atoms of the bond.
        default_radius (float): Defaults to 1. If ignore_outside is True, use this radius for atoms outside the bond.

    Returns:
        boolean array for each bond, whether it is occluded or not.
        If return_occluders is True, return the 2d boolean array marking occluding atoms instead.
    
    """
    periodic_table = Chem.GetPeriodicTable()

    if len(from_atoms) != len(to_atoms):
        raise ValueError("from_atoms and to_atoms must have the same length")
    if len(from_atoms) == 0:
        raise ValueError("from_atoms must not be empty")
    if len(potential_occluders) == 0:
        raise ValueError("potential_occluders must not be empty")
    

    # determine values required for the computation: positions and van-der-waals radii of the participating atoms
    if isinstance(from_atoms, np.ndarray) and (len(from_atoms.shape) == 2) and (from_atoms.shape[1] == 3):
        from_positions = from_atoms
        from_vdw = np.full(len(from_atoms), default_radius)
    else:
        from_positions = position(from_atoms)
        from_vdw = np.array([periodic_table.GetRvdw(atom.GetAtomicNum()) for atom in from_atoms])
    if isinstance(to_atoms, np.ndarray) and (len(to_atoms.shape) == 2) and (to_atoms.shape[1] == 3):
        to_positions = to_atoms
        to_vdw = np.full(len(to_atoms), default_radius)
    else:
        to_positions = position(to_atoms)
        to_vdw = np.array([periodic_table.GetRvdw(atom.GetAtomicNum()) for atom in to_atoms])
    if isinstance(potential_occluders, np.ndarray) and (len(potential_occluders.shape) == 2) and (potential_occluders.shape[1] == 3):
        occluders_positions = potential_occluders
        occluders_vdw = np.full(len(potential_occluders), default_radius)
    else:
        occluders_positions = position(potential_occluders)
        occluders_vdw = np.array([periodic_table.GetRvdw(atom.GetAtomicNum()) for atom in potential_occluders])

    # scale the van-der-waals radii by the vdw_factor
    occluders_vdw *= occluders_vdw_factor

    # for each of the bonds, determine the direction and the projection of the atoms onto the bond
    bond_directions = to_positions - from_positions; bond_directions /= np.linalg.norm(bond_directions, axis = 1)[:, None]

    from_projections = (bond_directions * from_positions).sum(axis = 1)
    to_projections = (bond_directions * to_positions).sum(axis = 1)
    occluders_projections = bond_directions @ occluders_positions.T


    # determine the caps of the projections a potential occluder might lie behind the bond but have a large enough vdw to still occlude
    from_projection_caps = from_projections + from_vdw
    to_projection_caps = to_projections - to_vdw
    if ignore_outside:
        occluders_filter = (occluders_projections > from_projection_caps[:, None]) & (occluders_projections < to_projection_caps[:, None])
    else:
        occluders_projections = np.clip(occluders_projections, from_projection_caps[:, None], to_projection_caps[:, None])
        occluders_filter = np.ones_like(occluders_projections, dtype = bool)

    # translate the projected points back to the original coordinate system and check if any vdw is large enough to occlude
    translated_projection = occluders_projections - from_projections[:, None]
    translated_projection = (translated_projection[:,:,None] @ bond_directions[:,None,:]) + from_positions[:,None,:]
    occluders_distance_to_projection = np.linalg.norm(translated_projection - occluders_positions[None,:,:], axis = -1)
    occluded = (occluders_distance_to_projection < occluders_vdw[None, :]) & occluders_filter
    if return_occluders:
        return occluded
    else:
        # if van der waals radii intersect, we assume no occlusion can happen
        return occluded.any(axis = -1) & (from_projection_caps < to_projection_caps)

def is_planar(atoms, threshold = 15):
    """
    Check if the given atoms are planar. The atoms are considered planar if the angle between any normal and an arbitrary reference normal is greater than threshold.

    Args:
        atoms (list): List of atoms to check for planarity.
        threshold (float): Defaults to 5. Threshold for the angle between the normal of the plane and the z-axis.

    Returns:
        bool: True if the atoms are planar, False otherwise.

    """
    atoms = position(atoms)
    atoms -= np.mean(atoms, axis = 0)
    reference_normal = np.cross(atoms[1] - atoms[0], atoms[2] - atoms[0])
    reference_normal /= np.linalg.norm(reference_normal)
    for i in range(3, len(atoms)):
        normal = np.cross(atoms[i] - atoms[0], atoms[0] - atoms[2])
        normal /= np.linalg.norm(normal)
        angle = np.arccos(np.clip(np.dot(normal, reference_normal), -1.0, 1.0)) * 180 / np.pi 
        if angle > threshold:
            return False
    return True


def score_overlap(l1, l2, occupation_max = 5, similarity_max = 5):
    """
    Scores overlap of l1 to l2, i.e. how much l1 lies within l2. Score is determined by:
        RMSD of:
            for each atom in l1:
                min(minimum distance to any atom in l2, occupation_max) + min(minimum distance to any atom of same element in l2, similarity_max)
        Devided by 2
    In order to only consider occupation or similarity, set the respective other max value to 0.

    Args:
        l1 (list): Rdkit Mol, list of atoms or atomkdtree of first molecule.
        l2 (list): Rdkit Mol, list of atoms or atomkdtree of second molecule.
        ocupation_max (float): Defaults to 5. Maximum distance to consider for occupation.
        similarity_max (float): Defaults to 5. Maximum distance to consider for similarity.

    Returns:
        float: Score of the overlap.

    """
    if isinstance(l1, AtomKDTree):
        l1 = l1
    elif isinstance(l1, list):
        l1 = AtomKDTree(l1)
    elif isinstance(l1, Chem.rdchem.Mol):
        l1 = AtomKDTree(l1.GetAtoms())
    else:
        raise TypeError("l1 must be a list of atoms, an rdkit molecule or an Atom")
    if isinstance(l2, AtomKDTree):
        l2 = l2
    elif isinstance(l2, list):
        l2 = AtomKDTree(l2)
    elif isinstance(l2, Chem.rdchem.Mol):
        l2 = AtomKDTree(l2.GetAtoms())
    else:
        raise TypeError("l2 must be a list of atoms, an rdkit molecule or an Atom")
    
    score = 0

    if not occupation_max == 0:
        # occupation check
        distance_matrix = l1.sparse_distance_matrix(l2, occupation_max).tocsr()
        distance_matrix.eliminate_zeros()
        min_distance = np.minimum.reduceat(distance_matrix.data, distance_matrix.indptr[:-1])
        min_distance[min_distance == 0] = occupation_max
        #occupations_scores = min_distance
        score += np.sqrt(np.mean(min_distance**2))

    # similarity check
    if not similarity_max == 0:
        min_distances = []
        for element in l1.trees.keys():
            if element not in l2.trees.keys():
                continue
            distance_matrix = l1.sparse_distance_matrix(l2, similarity_max, element).tocsr()
            distance_matrix.eliminate_zeros()
            min_distance = np.minimum.reduceat(distance_matrix.data, distance_matrix.indptr[:-1])
            min_distance[min_distance == 0] = occupation_max
            min_distances.extend(min_distance)
            #occupations_scores += min_distance
        min_distances = np.array(min_distances)
        score += np.sqrt(np.mean(min_distances**2))

    score = score / 2
    return score

def rotate_around_bond(mol, bond_start_atom_idx, bond_end_atom_idx, angle):
    """
    Rotates the atoms connected to the bond_end_atom by the given angle around the bond axis.
    Returns the modified molecule.

    Args:
        mol (RDKit.Mol): Molecule to rotate.
        bond_start_atom_idx (int): Index of the atom at the start of the bond.
        bond_end_atom_idx (int): Index of the atom at the end of the bond.
        angle (float): Angle to rotate the atoms by.
    
    Returns:
        RDKit.RWMol: Molecule with rotated atoms.
    """
    from .select import get_bond_connected_atoms
    from scipy.spatial.transform import Rotation as R
    # make mol RWMol
    # get bond and atoms
    bond_start_pos = position(mol.GetAtomWithIdx(bond_start_atom_idx))
    bond_end_pos = position(mol.GetAtomWithIdx(bond_end_atom_idx))
    bond_direction = bond_end_pos - bond_start_pos; bond_direction /= np.linalg.norm(bond_direction)
    rotation = R.from_rotvec(bond_direction * angle)
    # get atoms to rotate
    connected_atom_indices = get_bond_connected_atoms(mol, bond_end_atom_idx, bond_start_atom_idx)
    # rotate atoms
    original_position = position([mol.GetAtomWithIdx(int(i)) for i in connected_atom_indices])
    original_position -= bond_start_pos
    rotated_position = rotation.apply(original_position)
    rotated_position += bond_start_pos

    # apply positions to conformer
    mol = Chem.RWMol(mol)
    conf = next(mol.GetConformers())
    for atom_idx, pos in zip(connected_atom_indices, rotated_position):
        conf.SetAtomPosition(int(atom_idx), pos)
    return mol


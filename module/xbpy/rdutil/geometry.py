from rdkit.Chem import rdMolTransforms
from rdkit import Chem
import numpy as np
import logging
from scipy.spatial import cKDTree

class AtomKDTree():
    def __init__(self, atoms):
        vdw_radii = np.array([Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in atoms])
        self.trees = []
        for vdw_radius in np.unique(vdw_radii):
            tree = cKDTree(position(atoms)[vdw_radii == vdw_radius])
            self.trees.append((vdw_radius, tree))

    def query_ball_point(self, points, radius):
        indices = [set() for _ in range(len(points))]
        for vdw_radius, tree in self.trees:
            query_result = tree.query_ball_point(points, radius + vdw_radius)
            for i, result in enumerate(query_result):
                indices[i].update(result)
        return indices

def position(atom, conformer = 0):
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
        if not (type(conformer) == int):
            conformer = conformer
        else:
            conformer = atom[0].GetOwningMol().GetConformer(conformer)

        if len(atom) > 10: #TODO: 10 is currently a magic number, should be replaced by a more sensible value
            # assume its faster to get all positions at once and then slice
            
            molecular_positions = np.array(conformer.GetPositions())
            indices = [a.GetIdx() for a in atom]
            return molecular_positions[indices]
        else:
           return np.array([conformer.GetAtomPosition(a.GetIdx()) for a in atom])
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


def check_occlusion(from_atoms, to_atoms, potential_occluders, return_occluders = False, ignore_outside = True, default_radius = 1):
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


from rdkit.Chem import rdMolTransforms
from rdkit import Chem
import numpy as np


def position(atom, conformer = 0):
    """Return the position of the given atom.

    Args:
        atom (RDKit.Atom): Atom to get the position of. If list will return a list of positions.
        conformer (int): Defaults to 0. Conformer to get the position from.

    Returns:
        np.ndarray: 3xN Position of the given N atoms.

    """
    import numpy as np
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


def transform(mol, rotation, translation):
    transformation = np.eye(4)
    transformation[:3, :3] = rotation
    transformation[:3, 3] = translation
    rdMolTransforms.TransformConformer(mol.GetConformer(), transformation)


def check_occlusion(from_atoms, to_atoms, potential_occluders, return_occluders = False, ignore_outside = True):
    """
    Check if the bond between from_atoms and to_atoms is occluded by any of the atoms in potential_occluders.
    If return_occluders is True, return the occluding atoms. An atom is considered occluding the bond between two other atoms if its van-der-waals radius is larger than its distance to the bond.
    The bonds are assumed to be between the atoms at the same index in from_atoms and to_atoms.

    Args:
        from_atoms (list): List of atoms at the start of the bonds.
        to_atoms (list): List of atoms at the end of the bonds.
        potential_occluders (list): List of atoms that could occlude the bonds.
        return_occluders (bool): Defaults to False. If True, return a boolean array marking each occluding atom for each bond.
        ignore_outside (bool): Defaults to True. If True, only consider atoms that are between the two atoms of the bond.

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
    from_positions = position(from_atoms)
    from_vdw = np.array([periodic_table.GetRvdw(atom.GetAtomicNum()) for atom in from_atoms])
    to_positions = position(to_atoms)
    to_vdw = np.array([periodic_table.GetRvdw(atom.GetAtomicNum()) for atom in to_atoms])
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


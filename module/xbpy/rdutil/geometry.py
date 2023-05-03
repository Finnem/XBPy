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
        if len(atom) > 10: #TODO: 10 is currently a magic number, should be replaced by a more sensible value
            # assume its faster to get all positions at once and then slice
            molecular_positions = np.array(atom[0].GetOwningMol().GetConformer(conformer).GetPositions())
            indices = [a.GetIdx() for a in atom]
            return molecular_positions[indices]
        else:
            conformer = atom[0].GetOwningMol().GetConformer(conformer)
            return np.array([conformer.GetAtomPosition(a.GetIdx()) for a in atom])
    else:
        return np.array(atom.GetOwningMol().GetConformer(conformer).GetAtomPosition(atom.GetIdx()))
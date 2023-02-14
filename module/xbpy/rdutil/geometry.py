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
        molecular_positions = np.array(atom.GetOwningMol().GetConformer(conformer).GetPositions())
        indices = [a.GetIdx() for a in atom]
        return molecular_positions[indices]
    else:
        return np.array(atom.GetOwningMol().GetConformer(conformer).GetAtomPosition(atom.GetIdx()))
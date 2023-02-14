
def get_connected_atoms(atom, as_indices=False):
    """Return a list of atoms part of the connected component of the given atom.

    Args:
        atom (RDKit.Atom): Atom to get the connected atoms from.
        as_indices (bool): Defaults to False. If True, return the indices of the atoms instead of the atoms themselves.

    Returns:
        list(RDKit.Atom): List of atoms part of the connected component of the given atom.

    """

    molecule = atom.GetOwningMol()
    # we keep a list and a set to keep track of the order of when we saw the atoms but have fast lookup
    connected_atom_ids = [atom.GetIdx()]
    seen = set(connected_atom_ids) 
    to_check = list(connected_atom_ids)
    while len(to_check) > 0:
        atom_idx = to_check.pop(0)
        for neighbor in molecule.GetAtomWithIdx(atom_idx).GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in seen:
                seen.add(neighbor_idx)
                to_check.append(neighbor_idx)
                connected_atom_ids.append(neighbor_idx)
    if as_indices:
        return connected_atom_ids
    else:
        return [molecule.GetAtomWithIdx(atom_idx) for atom_idx in connected_atom_ids]

    
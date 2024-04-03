from rdkit.Chem import GetPeriodicTable
from rdkit.Chem.rdchem import Atom
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



def get_small_molecules(structure, min_atoms = 6, max_atoms = 200):
    """
        By default returns all small molecules with 6 to 200 atoms.
    """

    relevant_atom_ids = set()
    for atom in structure.GetAtoms():
        if atom.GetAtomicNum() != 1:
            relevant_atom_ids.add(atom.GetIdx())

    small_molecule_indices = []
    while len(relevant_atom_ids) > 0:
        atom_idx = relevant_atom_ids.pop()
        connected_atom_indices = get_connected_atoms(structure.GetAtomWithIdx(atom_idx), as_indices=True)
        if (len(connected_atom_indices) >= min_atoms) and (len(connected_atom_indices) <= max_atoms):
            small_molecule_indices.append(connected_atom_indices)
        relevant_atom_ids.difference_update(connected_atom_indices)
    
    return [[structure.GetAtomWithIdx(a_idx) for a_idx in (atom_indices)] for atom_indices in small_molecule_indices]

def vdw_radius(atom):
    """Return the van der Waals radius of the given atom.

    Args:
        atom (RDKit.Atom): Atom to get the van der Waals radius from. Can be a list of atoms.

    Returns:
        float: van der Waals radius of the given atom.

    """
    try:
        atom = list(atom)
    except TypeError:
        atom = [atom]

    atomic_number = [a.GetAtomicNum() if type(a) == Atom else GetPeriodicTable().GetAtomicNumber(a) for a in atom]

    return [GetPeriodicTable().GetRvdw(a) for a in atom]
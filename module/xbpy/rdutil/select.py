from rdkit.Chem import GetPeriodicTable
from rdkit.Chem.rdchem import Atom
from rdkit import Chem


def bond_independent_substructure_matches(mol, query):
    mol_copy = Chem.Mol(mol)
    query_copy = Chem.Mol(query)
    for bond in mol_copy.GetBonds():
        bond.SetBondType(Chem.BondType.UNSPECIFIED)
    for bond in query_copy.GetBonds():
        bond.SetBondType(Chem.BondType.UNSPECIFIED)
    return mol_copy.GetSubstructMatches(query_copy)

def select_atom(possible_atoms, neighborhood = None, element = None, return_first = True):
    """
    Returns a atom (if return first) or all atoms that have the given amount of neighbors.

    Args:
        possible_atoms (RDKit.Mol or iterable of atoms): List of Atoms or Molecule to get the atoms from.
        neighborhood (dict): Dictionary mapping atom symbols to their number of occurences that is beeing sought.
        return_first (bool): Defaults to True. If True, return the first atom found. If False, return all atoms found.

    Returns:
        list(RDKit.Atom) or RDKit.Atom: List of atoms that have the given amount of neighbors.

    """
    atoms = []
    if neighborhood is None:
        neighborhood = {}
    if hasattr(possible_atoms, "GetAtoms"):
        possible_atoms = possible_atoms.GetAtoms()
    for atom in possible_atoms:
        invalid = False
        if not element is None:
            if atom.GetSymbol() != element:
                continue
        for symbol, count in neighborhood.items():
            for n in atom.GetNeighbors():
                if n.GetSymbol() == symbol:
                    count -= 1
            if count != 0:
                invalid = True
                break
        if not invalid:
            if return_first:
                return atom
            atoms.append(atom) 
    return atoms

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
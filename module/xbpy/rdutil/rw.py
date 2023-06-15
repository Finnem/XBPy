import numpy as np


def remove_atoms(mol, atoms):
    """Remove the given atoms from the molecule.
    
    Args:
        mol (RDKit.Mol): Molecule to remove the atoms from.
        atoms (list(RDKit.Atom or int)): List of atoms to remove. If list of ints, the atoms with the given indices will be removed.
            Passing indices will generally be faster.        
    Returns:
        RDKit.Mol: Molecule with the given atoms removed.
        
    """
    import rdkit.Chem as Chem
    from .geometry import position

    indices = []
    for atom in atoms:
        if isinstance(atom, int):
            indices.append((atom, None))
        else:
            indices.append((atom.GetIdx(), atom))
    # we need to remove the atoms in reverse order to not mess up the indices
    indices.sort(reverse=True)
    mol = Chem.RWMol(mol)
    removed = set()
    for atom_idx, atom in indices:
        if atom_idx in removed:
            continue
        if atom is not None:
            other_atom = mol.GetAtomWithIdx(atom_idx)
            if (other_atom.GetSymbol() != atom.GetSymbol()) or (not np.allclose(position(other_atom), position(atom))):
                raise ValueError("The atom at index {} is not the same as the given atom.".format(atom_idx))
            mol.RemoveAtom(atom_idx)
        else:
            mol.RemoveAtom(atom_idx)
    return mol.GetMol()

def keep_atoms(mol, atoms):
    """Remove all atoms except the given atoms from the molecule.
    
    Args:
        mol (RDKit.Mol): Molecule to remove the atoms from.
        atoms (list(RDKit.Atom or int)): List of atoms to keep. If list of ints, the atoms with the given indices will be kept.
            Passing indices will generally be faster.
            
    Returns:
        RDKit.Mol: Molecule with all atoms except the given atoms removed.
        
    """
    all_atom_indices = set([atom.GetIdx() for atom in mol.GetAtoms()])
    indices = set()
    for atom in atoms:
        if isinstance(atom, int):
            indices.add(atom)
        else:
            indices.add(atom.GetIdx())
    to_remove = all_atom_indices - indices
    return remove_atoms(mol, to_remove)


def copy_props(mol_from, mol_to, replace_dict = None):
    if replace_dict is None:
        replace_dict = {}

    for key, prop in mol_from.GetPropsAsDict().items():
        if key in replace_dict:
            key = replace_dict[key]
        if np.issubdtype(type(prop), np.integer):
            mol_to.SetIntProp(key, prop)
        elif np.issubdtype(type(prop), np.floating):
            mol_to.SetDoubleProp(key, prop)
        elif isinstance(type(prop), str):
            mol_to.SetProp(key, prop)

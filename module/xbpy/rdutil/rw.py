import numpy as np
from rdkit.Chem import AllChem
import logging
from rdkit.Geometry import Point3D
from rdkit import Chem
from ..math.geometry import calculate_angle
from .geometry import position
from itertools import combinations
from rdkit.Chem import rdmolops
from .bond_order_inference import correct_bond_orders as correct_bond_orders_func
from rdkit.Chem.PropertyMol import PropertyMol

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
        if np.issubdtype(type(atom), np.integer):
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
        if np.issubdtype(type(atom), np.integer):
            indices.add(atom)
        else:
            indices.add(atom.GetIdx())
    to_remove = all_atom_indices - indices
    return remove_atoms(mol, to_remove)

ideal_tetraheder = np.array([[0,0,0], [0.0, 0.0, 1.0], [0.0, 0.8944271909999159, -0.4472135954999579], [-0.816496580927726, -0.408248290463863, -0.408248290463863], [0.816496580927726, -0.408248290463863, -0.408248290463863]])

def add_explicit_hydrogens(mol):
    """Crude function to add hydrogens to a molecule in a simple fashion.
    Currently only works for the organic elements.
    
    Args:
        mol (RDKit.Mol): Molecule to add hydrogens to.
        
    Returns:
        RDKit.Mol: Molecule with explicit hydrogens added.
        
    """
    from rdkit import Chem
    from .geometry import position
    from ..math.geometry import rigid_transform, apply_transform
    from rdkit.Chem.rdmolops import RemoveHs

    mol = RemoveHs(mol, implicitOnly=True)
    mol = Chem.RWMol(mol)
    
    periodic_table = Chem.GetPeriodicTable()
    for atom in list(mol.GetAtoms()):
        possible_valence = periodic_table.GetDefaultValence(atom.GetAtomicNum())

        neighbors = list(atom.GetNeighbors())
        cur_valence = atom.GetExplicitValence()
        missing_hydrogens =  possible_valence - cur_valence

        if missing_hydrogens > 0:
            neighbor_directions = np.array([position(neighbor) - position(atom) for neighbor in neighbors]); neighbor_directions /= np.linalg.norm(neighbor_directions, axis=1)[:,None]
            neighbor_directions = np.concatenate([np.zeros((1, 3)), neighbor_directions])
            
            hydrogen_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() == 1]
            hydrogen_distance = 1.0 if not atom.GetSymbol() == "C" else 1.1
            if len(hydrogen_neighbors) > 0:
                hydrogen_distance = np.mean([np.linalg.norm(position(neighbor) - position(atom)) for neighbor in hydrogen_neighbors])
            if possible_valence == 4:
                # sp3 hybridization => match tetrahedral geometry to neighbors

                    
                    R, t  = rigid_transform(ideal_tetraheder[:len(neighbor_directions)], neighbor_directions)
                    new_positions = apply_transform(ideal_tetraheder[len(neighbor_directions):], R, t)
                    for new_position in new_positions:
                        new_position = new_position * hydrogen_distance + position(atom)
                        mol.AddAtom(Chem.Atom(1))
                        new_idx = mol.GetNumAtoms() - 1
                        mol.AddBond(atom.GetIdx(), new_idx, Chem.BondType.SINGLE)
                        mol.GetConformer().SetAtomPosition(new_idx, new_position)

    return mol
                  


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
        else:
            mol_to.SetProp(str(key), prop)



def build_molecule(positions, elements, bond_indices, bond_orders, formal_charges=None, hybridizations=None, sanitize=True):
    """
    Constructs an RDKit molecule from atom positions, element symbols, bond indices, bond orders,
    and optionally formal charges and hybridization states.

    Args:
        positions (numpy.ndarray): A NumPy array of shape (N, 3) representing the (x, y, z) coordinates of N atoms.
        elements (list): List of element symbols (e.g., ['C', 'O', 'N']) corresponding to the atoms.
        bond_indices (list): List of tuples representing the bonded atom indices. Bonds are assumed to be between atoms at the same index in `from_atoms` and `to_atoms`.
        bond_orders (list): List of bond orders as strings (e.g., "SINGLE", "DOUBLE", "TRIPLE", "AROMATIC", "ONEANDAHALF").
        formal_charges (list, optional): List of formal charges corresponding to each atom. If not provided, formal charges default to 0 for all atoms.
        hybridizations (list, optional): List of hybridization states for each atom as strings (e.g., "SP", "SP2", "SP3"). If not provided, hybridization defaults to unspecified for all atoms.
        sanitize (bool, optional): Whether to sanitize the molecule after construction. Defaults to True.

    Returns:
        rdkit.Chem.rdchem.Mol: An RDKit molecule object with the specified atoms, bonds, and coordinates.
    
    Raises:
        ValueError: If the positions array does not have the shape (N, 3).
    
    Example:
        positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        elements = ['C', 'O', 'N']
        bond_indices = [(0, 1), (1, 2)]
        bond_orders = ["DOUBLE", "SINGLE"]
        formal_charges = [0, -1, 0]
        hybridizations = ["SP2", "SP3", "SP2"]

        molecule = build_molecule(positions, elements, bond_indices, bond_orders, formal_charges, hybridizations)
    """
    # Ensure positions is a NumPy array
    if not isinstance(positions, np.ndarray):
        positions = np.array(positions)
    
    # Ensure the array has the right shape
    if positions.shape[1] != 3:
        raise ValueError("Positions array must have shape (N, 3) for 3D coordinates")
    
    # Create an editable molecule
    mol = Chem.RWMol()
    
    # Hybridization state mapping
    hybridization_mapping = {
        "SP": Chem.HybridizationType.SP,
        "SP2": Chem.HybridizationType.SP2,
        "SP3": Chem.HybridizationType.SP3,
        "SP3D": Chem.HybridizationType.SP3D,
        "SP3D2": Chem.HybridizationType.SP3D2
    }
    
    # Add atoms to the molecule and assign optional formal charges and hybridization states
    for i, elem in enumerate(elements):
        atom = Chem.Atom(elem)
        
        # Set formal charge if provided
        if formal_charges is not None:
            atom.SetFormalCharge(int(formal_charges[i]))
        
        # Set hybridization if provided
        if hybridizations is not None:
            hybrid_state = hybridization_mapping.get(hybridizations[i], Chem.HybridizationType.UNSPECIFIED)
            atom.SetHybridization(hybrid_state)
        
        mol.AddAtom(atom)
    
    # Map bond order strings to RDKit bond types
    bond_order_mapping = {
        "SINGLE": Chem.BondType.SINGLE,
        "DOUBLE": Chem.BondType.DOUBLE,
        "TRIPLE": Chem.BondType.TRIPLE,
        "AROMATIC": Chem.BondType.AROMATIC,
        "ONEANDAHALF": Chem.BondType.AROMATIC  # Map ONEANDAHALF to AROMATIC bond type
    }
    
    # Add bonds to the molecule
    for (start_idx, end_idx), bond_order in zip(bond_indices, bond_orders):
        bond_type = bond_order_mapping.get(bond_order, Chem.BondType.SINGLE)
        mol.AddBond(int(start_idx), int(end_idx), bond_type)
    
    # Create a conformer and set atom positions
    conf = Chem.Conformer(len(positions))
    for i, pos in enumerate(positions):
        x, y, z = pos
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf)
    
    # Finalize the molecule
    mol.UpdatePropertyCache(strict=False)
    if sanitize:
        Chem.SanitizeMol(mol)
    
    return mol


def kekulize_mol(mol):
    from rdkit.Chem import rdmolops
    from rdkit.Chem import rdchem
    rdmolops.FastFindRings(mol)

    mol = mol.GetMol()
    # remove implicit hydrogens
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    atom_to_ring = {atom_idx : set() for atom_idx in range(mol.GetNumAtoms())}
    for ring in atom_rings:
        for atom_idx in ring:
            atom_to_ring[atom_idx].add(ring)
    possibly_aromatic = set()
    for atom in mol.GetAtoms():
        if atom.GetNumImplicitHs() > 0:
            possibly_aromatic.add(atom.GetIdx())
        atom.SetNoImplicit(True)  # Disable automatic calculation of implicit hydrogens
    for bond in mol.GetBonds():
        if (bond.GetBondType() == rdchem.BondType.SINGLE) and \
        ((bond.GetBeginAtomIdx() in possibly_aromatic) or (bond.GetEndAtomIdx() in possibly_aromatic)) and \
        (atom_to_ring[bond.GetBeginAtomIdx()] & atom_to_ring[bond.GetEndAtomIdx()]):
            bond.SetBondType(rdchem.BondType.AROMATIC)

    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # if kekulize didnt catch the ring, its not aromatic
    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.AROMATIC:
            if not bond.GetIsAromatic():
                bond.SetBondType(rdchem.BondType.SINGLE)

    mol.UpdatePropertyCache(strict=True)
    return mol

def proximity_bond(mol, correct_bond_orders = True, allow_hydrogen_gas = False, try_full = True, as_property_mol = False):
    """
    Infer bonds from the proximity of atoms. This is done by checking the distance of each atom to each other atom.
    If the distance is smaller than the sum of the van-der-waals radii, a bond is inferred.
    This method is not transitive, i.e. if atom A is close to atom B and atom B is close to atom C, A and C are not necessarily close.
    A new molecule is created and returned.

    Args:
        mol (RDKit.Mol): Molecule to infer bonds for.
        kekulize (bool): Defaults to True. If True, molecule will be kekulized after adding bonds.
        handle_valency (str): Defaults to "auto". If "auto", valency violations will be resolved by removing the most unlikely bond. If "ignore", valency violations will be ignored. If "raise", valency violations will raise an error.
        allow_hydrogen_gas (bool): Defaults to False. If True, hydrogen gas will be allowed. This means that hydrogen atoms will not be considered to have valency violations.

    Returns:
        RDKit.Mol: Molecule with inferred bonds.

    """
    from rdkit.Chem import rdchem
    from rdkit.Chem import rdMolTransforms
    from rdkit.Chem import rdmolops
    from .util import bond_lengths, ideal_bond_lengths
    from .geometry import position
    from scipy.spatial import cKDTree

    # get the positions of the atoms
    positions = position(mol)
    periodic_table = Chem.GetPeriodicTable()
    elements = np.array([atom.GetSymbol() for atom in mol.GetAtoms()])
    # create a kdtree for the atoms
    tree = cKDTree(positions)

    # query the tree for atoms that are closer than maximum sensible bond length
    d_matrix = tree.sparse_distance_matrix(tree, max_distance = 2.5, output_type = 'coo_matrix')

    # get bond_distances for each pair of atoms
    bond_distances = np.array([bond_lengths[elements[i]][elements[j]] for i, j in zip(d_matrix.row, d_matrix.col)])
    
    cond = (d_matrix.data < bond_distances[:,0] * 1.05) & (d_matrix.row < d_matrix.col)
    close_atoms = np.vstack([d_matrix.row[cond], d_matrix.col[cond]]).T
    mol = Chem.RWMol(mol)
    # create a bond for each pair of close atoms
    for idx, (i, j) in enumerate(close_atoms):
        if not allow_hydrogen_gas:
            if (elements[i] == "H") and (elements[j] == "H"):
                continue
        dist = d_matrix.data[cond][idx]
        if dist < bond_distances[cond][idx][3]:
            if mol.GetBondBetweenAtoms(int(i), int(j)) is None:
                mol.AddBond(int(i), int(j), rdchem.BondType.TRIPLE)
            else:
                mol.GetBondBetweenAtoms(int(i), int(j)).SetBondType(rdchem.BondType.TRIPLE)
        elif dist < bond_distances[cond][idx][2]:
            if mol.GetBondBetweenAtoms(int(i), int(j)) is None:
                mol.AddBond(int(i), int(j), rdchem.BondType.DOUBLE)
            else:
                mol.GetBondBetweenAtoms(int(i), int(j)).SetBondType(rdchem.BondType.DOUBLE)
        elif dist < bond_distances[cond][idx][1]:
            bond_order = rdchem.BondType.SINGLE#ONEANDAHALF
            if mol.GetBondBetweenAtoms(int(i), int(j)) is None:
                mol.AddBond(int(i), int(j), bond_order)
            else:
                mol.GetBondBetweenAtoms(int(i), int(j)).SetBondType(bond_order)
        else:
            if mol.GetBondBetweenAtoms(int(i), int(j)) is None:
                mol.AddBond(int(i), int(j), rdchem.BondType.SINGLE)
            else:
                mol.GetBondBetweenAtoms(int(i), int(j)).SetBondType(rdchem.BondType.SINGLE)

    # check for valency violations
    mol.UpdatePropertyCache(strict=False) 
    if correct_bond_orders:
        mol = correct_bond_orders_func(mol, try_full=try_full)

    return PropertyMol(mol) if as_property_mol else mol.GetMol()
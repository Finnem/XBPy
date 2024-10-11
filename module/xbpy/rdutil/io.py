
import logging
from tqdm import tqdm
import numpy as np
import re
from rdkit import Chem
from rdkit.Chem.rdmolops import RenumberAtoms
from .util import jump_to_nth_last_line


def determine_molecule_paths(paths, recursive = True):
    import os 
    import glob
    molecule_paths = []
    for path in paths:
        remaining_paths = set([path])
        seen = set()
        path_molecules = []
        while remaining_paths:
            p = remaining_paths.pop()
            for f in glob.glob(p):
                if os.path.isdir(f):
                    # If recursive we do a random walk through the directory tree
                    if (f not in seen) and recursive:
                        remaining_paths.add(f)
                else:
                    molecule_paths.append(f)
                    
            seen.add(p)
        if len(molecule_paths) == 0:
            logging.warning("No molecules found at path(s) {}. Make sure the file(s) exists.".format(path))
    return molecule_paths

def read_molecules(path, recursive = True, store_path = False, reference_molecule = None, removeHs=False, sanitize=False, proximityBonding=True, reset_index = False, *args, **kwargs):
    """Read molecules as RDK molecules from a single pdb, mol or sdf file, or from a directory /multiple directories of such files.
        Coordinates from an xyz file are converted to a RDKit molecule using the reference molecule as a template.
    
    Args:
        path (str or list(str)): Path to a file or directory. May be a list of paths. A path may also describe a glob pattern.
        recursive (bool): Defaults to True. If True, recursively search directories for files. 
        *args: Additional arguments to pass to the RDKit reader.
        **kwargs: Additional keyword arguments to pass to the RDKit reader.

    Yields:
        RDKit.Mol: RDKit molecules.

    """

    import os
    import glob

    if isinstance(path, str):
        path = [path]
    paths = path

    molecules = []
    # Loop over all paths, seperated from tree walk to allow for error detection (emtpy paths)
    molecule_paths = determine_molecule_paths(paths, recursive=recursive)
    for path in sorted(molecule_paths):
        path_molecules = _read_molecules_file(path, store_path, reference_molecule, removeHs=removeHs, sanitize=sanitize, proximityBonding=proximityBonding, *args, **kwargs)
        if len(path_molecules) == 0:
            logging.warning("No molecules found at path(s) {}. Make sure the file(s) exists.".format(path))
        else:
            for mol in path_molecules:
                if reset_index:
                    from ..morgan import unique_index
                    new_indices = unique_index(mol)
                    ordering = np.argsort(new_indices)
                    mol = RenumberAtoms(mol, [int(i) for i in ordering])
                yield mol


def get_num_molecules(paths, recursive = True, *args, **kwargs):
    """
    Determine Number of molecules in a file or directory. Number is determined as lazily as possible.

    Args:
        path (str or list(str)): Path to a file or directory. May be a list of paths. A path may also describe a glob pattern.
        recursive (bool): Defaults to True. If True, recursively search directories for files.

    Returns:
        int: Number of molecules.

    """
    if isinstance(paths, str):
        paths = [paths]
    molecule_paths = determine_molecule_paths(paths, recursive=recursive)

    num_molecules = 0
    for paths in molecule_paths:
        num_molecules += _get_num_mols(paths, *args, **kwargs)

    return num_molecules


def _get_num_mols(path, *args, **kwargs):
    """Try to determine number of molecules for given path as fast as possible.

    Args:
        path (str): Path to a file or directory.

    Returns:
        int: Number of molecules.
    """

    import os
    if os.path.splitext(path)[1] == ".sdf": # for sdf we can jump to the end of the file and read the numbering of an attribute
        attribute_line = jump_to_nth_last_line(path, 4)
        match = re.search(r'\((.*?)\)', attribute_line)
        if match:
            return int(match.group(1))
        else: # no attributes present, then we have to count the $$$$ occurences
            with open(path, "r") as f:
                return sum(1 for line in f if line[:4] == "$$$$")
    else:
        return 1




def _read_molecules_file(path, store_path = True, reference_molecule = None, proximityBonding=True, *args, **kwargs):
    """Read molecules as RDK molecules from a single pdb, mol or sdf file. A xyz file is converted to a RDKit molecule using the reference molecule as a template.

    Args:
        path (str): Path to a file or directory. May be a list of paths. A path may also describe a glob pattern.

    Returns:
        list(RDKit.Mol): List of RDKit molecules.

    """

    import os
    import rdkit.Chem as Chem
    from rdkit.Chem import rdmolfiles

    molecules = []
    if os.path.splitext(path)[1] == ".pdb":
        if not "flavor" in kwargs:
            kwargs["flavor"] = 4
        molecules = [Chem.rdmolfiles.MolFromPDBFile(path, proximityBonding=proximityBonding, *args, **kwargs)]
    elif os.path.splitext(path)[1] == ".mol2":
        molecules = [Chem.rdmolfiles.MolFromMol2File(path, *args, **kwargs)]
    elif os.path.splitext(path)[1] == ".mol":
        molecules = [Chem.rdmolfiles.MolFromMolFile(path, *args, **kwargs)]
    elif os.path.splitext(path)[1] == ".sdf":
        molecules = Chem.rdmolfiles.SDMolSupplier(path, *args, **kwargs)
    elif os.path.splitext(path)[1] == ".xyz":
        molecules = read_molecules_xyz(path, reference_molecule, proximityBonding=proximityBonding, *args, **kwargs)
    else:
        try:
            molecules = read_coord_file(path, reference_molecule, proximityBonding=proximityBonding, *args, **kwargs)
        except KeyboardInterrupt:
            raise
        except:
            logging.warning(f"File format not recognized for {path}.")

    if store_path:
        new_mols = []
        for molecule in molecules:
            prop_names = molecule.GetPropNames()
            i = 1
            if "path" in prop_names:
                while f"path_{i}" in prop_names:
                    i += 1
                molecule.SetProp(f"path_{i}", path)
            else:
                molecule.SetProp("path", path)
            new_mols.append(molecule)
        molecules = new_mols

    if reference_molecule is not None:
        new_mols = []
        for molecule in molecules:
            new_mol = Chem.Mol(reference_molecule)
            for prop in reference_molecule.GetPropNames():
                new_mol.ClearProp(prop)
            for prop in molecule.GetPropNames():
                new_mol.SetProp(prop, molecule.GetProp(prop))
            for atom in new_mol.GetAtoms():
                conformer = new_mol.GetConformer()
                xyz_conformer = molecule.GetConformer()
                conformer.SetAtomPosition(atom.GetIdx(), xyz_conformer.GetAtomPosition(atom.GetIdx()))
                atom.SetAtomicNum(molecule.GetAtomWithIdx(atom.GetIdx()).GetAtomicNum())
            new_mols.append(new_mol)
        molecules = new_mols

    return molecules


def read_molecules_xyz(path, reference_molecule =None , proximityBonding = True, *args, **kwargs):
    """Read molecules as RDK molecules from a single xyz file. Coordinates are converted to a RDKit molecule using the reference molecule as a template.

    Args:
        path (str): Path to a file or directory. May be a list of paths. A path may also describe a glob pattern.
        reference_molecule (RDKit.Mol): Reference molecule to use as a template.

    Returns:
        list(RDKit.Mol): List of RDKit molecules.

    """

    import os
    import rdkit.Chem as Chem

    if os.path.splitext(path)[1] == ".xyz":
        xyz_molecule = Chem.rdmolfiles.MolFromXYZFile(path)
        # we infer the bonds if requested
        if proximityBonding:
            mol = proximity_bond(mol)
            rdmolops.RemoveHs(mol, implicitOnly=False, updateExplicitCount=False, sanitize=False)

        # finally we use the reference molecule to copy the bonds if sensible
        if not (reference_molecule is None):
            if not proximityBonding:
                raise ValueError("Reference molecule provided but proximity bonding not requested.")
            if not type(reference_molecule) == list:
                reference_molecule = [reference_molecule]
            mol = reduce_to_templates(mol, reference_molecule)

    return [mol]


def reduce_to_templates(mol, reference_molecules, copy_bond_order = True):
    """ 
    Searches reference molecules as substructures in the input molecule and copies the bonds correspondingly.
    Removes any atoms that are not part of the substructure.
    Args:
        mol (RDKit.Mol): RDKit molecule.
        reference_molecules (list(RDKit.Mol)): List of RDKit molecules to use as templates.
        copy_bond_order (bool): Defaults to True. If True, the bond order is copied as well.

    Returns:
        RDKit.Mol: RDKit molecule with bonds copied from the reference molecules.
    """

    from rdkit.Chem import AllChem

    mol = Chem.RWMol(mol)
    bonds_to_remove = set()

    mol_copy = Chem.Mol(mol)

    # prepare bond orders to be matched
    for bond in mol_copy.GetBonds():
        bonds_to_remove.add((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        bond.SetBondType(Chem.BondType.SINGLE)

    reverse_map = {}
    mapped = set()

    # determine substructure matches
    for reference_molecule in reference_molecules:
        ref_copy = Chem.Mol(reference_molecule)
        for bond in ref_copy.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
        for match in mol.GetSubstructMatches(ref_copy):
            mapped.update(match)
            for i, bond in enumerate(ref_copy.GetBonds()):
                to_keep = (match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()])
                reverse_map[to_keep] = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                if to_keep[1] < to_keep[0]:
                    to_keep = (to_keep[1], to_keep[0])
                bonds_to_remove.difference_update([to_keep])

    # remove any bond that was not found in the reference molecules
    for bond in bonds_to_remove:
        mol.RemoveBond(*bond)

    # assign bond orders from reference molecule
    if copy_bond_order:
        for bond in mol.GetBonds():
            ref_bond = reverse_map.get((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            if ref_bond is not None:
                ref_bond = ref_copy.GetBondBetweenAtoms(ref_bond[0], ref_bond[1])
                bond.SetBondType(ref_bond.GetBondType())

    # remove unmapped_atoms
    unmapped_atoms = set(range(mol.GetNumAtoms())) - mapped
    unmapped_atoms = sorted(list(unmapped_atoms), reverse=True)
    for atom in unmapped_atoms:
        mol.RemoveAtom(atom)

    return mol.GetMol()



def read_coord_file(path, reference_molecule = None, proximityBonding = True, to_angstrom = True, *args, **kwargs):

    """
    Method to read molecules from a turbomole coord file. User-defined bonds are not supported.
    Works by first converting the coord file to a xyz block, infer the bonds and potentially cull bonds to the reference molecule.

    Args:
        path (str): Path to the coord file.
        reference_molecule (RDKit.Mol): Reference molecule to use as a template.
        proximityBonding (bool): Defaults to False. If True, infer bonds from the xyz block. If False, use the reference molecule as a template.
        to_angstrom (bool): Defaults to True. If True, the coordinates are in bohr units and will be converted to angstrom.
        *args: Additional arguments to pass to the RDKit reader.
        **kwargs: Additional keyword arguments to pass to the RDKit reader.

    Returns:
        RDKit.Mol: List of RDKit molecules.
    """
    from .geometry import proximity_bond
    from rdkit.Chem import rdmolops


    # first we read the coord file using standard python IO
    lines = []

    # coord files are usually in bohr units, so we need to scale them to angstrom
    if to_angstrom:
        scaling_factor = 0.529177
    else:  
        scaling_factor = 1

    # we basically parse out the coords and elements, xyz only differ by placing the element first.
    coord_section = False
    for line in open(path, "r"):
        if "$coord" in line:
            coord_section = True
        elif coord_section:
            if line.startswith("$"):
                break
            else:
                lineparts = line.split()
                element = lineparts[3]
                element = element[0].upper() + element[1:].lower()
                other_lineparts = '\t'.join([f"{float(c) * scaling_factor}" for c in lineparts[:3]])
                new_line = f"{element}\t{other_lineparts}\n"
                lines.append(new_line)

    # finally we insert the number of atoms as based on the definition
    lines.insert(0, f"{len(lines)}\n\n")
    xyz_block = "".join(lines)

    # now we can construct the molecule by passing in the constructed xyz block to rdkit
    mol = Chem.MolFromXYZBlock(xyz_block)

    # and we infer the bonds if requested
    if proximityBonding:
        mol = proximity_bond(mol)
        rdmolops.RemoveHs(mol, implicitOnly=False, updateExplicitCount=False, sanitize=False)

    # finally we use the reference molecule to copy the bonds if sensible
    if not (reference_molecule is None):
        if not proximityBonding:
            raise ValueError("Reference molecule provided but proximity bonding not requested.")
        if not type(reference_molecule) == list:
            reference_molecule = [reference_molecule]
        mol = reduce_to_templates(mol, reference_molecule)

    return [mol]





def write_as_batches(molecules, batch_size, type = "xyz"):
    ...

import logging
from tqdm import tqdm
import numpy as np
from rdkit.Chem.rdmolops import RenumberAtoms

def read_molecules(path, recursive = True, store_path = False, reference_molecule = None, removeHs=False, sanitize=False, proximityBonding=False, reset_index = False, *args, **kwargs):
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


def _read_molecules_file(path, store_path = True, reference_molecule = None, proximityBonding=False, *args, **kwargs):
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
        raise ValueError("File format not recognized.")

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


def read_molecules_xyz(path, reference_molecule, proximityBonding = False, *args, **kwargs):
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
        if reference_molecule is None:
            #logging.warning("No reference molecule provided. Inferring molecule from xyz only.")
            if proximityBonding:
                molecule = Chem.MolFromPDBBlock(Chem.MolToPDBBlock(xyz_molecule), proximityBonding=True, sanitize=False, removeHs=False)
                logging.warning("Inferring bonds from xyz file. This may lead to incorrect results.")
            else:
                molecule = xyz_molecule
            try:
                molecule.UpdatePropertyCache()
                Chem.rdmolops.AssignAtomChiralTagsFromStructure(molecule)
                Chem.rdmolops.AssignStereochemistry(molecule)
            except Chem.AtomValenceException:
                if not kwargs.get("sanitize", False):
                    logging.warning("Could not infer stereochemistry from xyz file. This may be due to the molecule not being sanitized.")
                else:
                    raise

        else:
            molecule = Chem.Mol(reference_molecule)
            for prop in reference_molecule.GetPropNames():
                molecule.ClearProp(prop)
            for atom in molecule.GetAtoms():
                conformer = molecule.GetConformer()
                xyz_conformer = xyz_molecule.GetConformer()
                conformer.SetAtomPosition(atom.GetIdx(), xyz_conformer.GetAtomPosition(atom.GetIdx()))
                atom.SetAtomicNum(xyz_molecule.GetAtomWithIdx(atom.GetIdx()).GetAtomicNum())


    else:
        raise ValueError("File format not recognized.")
    


    return [molecule]


def write_as_batches(molecules, batch_size, type = "xyz"):
    ...

import logging

def read_molecules(path, recursive = True, store_path = False, reference_molecule = None, maximum = None, *args, **kwargs):
    """Read molecules as RDK molecules from a single pdb, mol or sdf file, or from a directory /multiple directories of such files.
        Coordinates from an xyz file are converted to a RDKit molecule using the reference molecule as a template.
    
    Args:
        path (str or list(str)): Path to a file or directory. May be a list of paths. A path may also describe a glob pattern.
        recursive (bool): Defaults to True. If True, recursively search directories for files. 
        *args: Additional arguments to pass to the RDKit reader.
        **kwargs: Additional keyword arguments to pass to the RDKit reader.

    Returns:
        list(RDKit.Mol): List of RDKit molecules.

    """

    import os
    import glob

    if isinstance(path, str):
        path = [path]
    paths = path

    molecules = []
    # Loop over all paths, seperated from tree walk to allow for error detection (emtpy paths)
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
                    used_maximum = None if maximum is None else maximum - (len(molecules) + len(path_molecules))
                    path_molecules.extend(_read_molecules_file(f, store_path, reference_molecule, maximum = used_maximum, *args, **kwargs))
                    
            seen.add(p)

        if len(path_molecules) == 0:
            logging.warning("No molecules found at path(s) {}. Make sure the file(s) exists.".format(path))
        else:
            molecules.extend(path_molecules)
    return molecules


def _read_molecules_file(path, store_path = True, reference_molecule = None, maximum = None, *args, **kwargs):
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
        molecules = [Chem.rdmolfiles.MolFromPDBFile(path, *args, **kwargs)]
    elif os.path.splitext(path)[1] == ".mol":
        molecules = [Chem.rdmolfiles.MolFromMolFile(path, *args, **kwargs)]
    elif os.path.splitext(path)[1] == ".sdf":
        molecules = Chem.rdmolfiles.SDMolSupplier(path, *args, **kwargs)
    elif os.path.splitext(path)[1] == ".xyz":
        molecules = read_molecules_xyz(path, reference_molecule, *args, **kwargs)
    else:
        raise ValueError("File format not recognized. RDKit does not properly support mol2 files.")

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


def read_molecules_xyz(path, reference_molecule, *args, **kwargs):
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
            molecule = xyz_molecule
            molecule.UpdatePropertyCache()
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(molecule)
            Chem.rdmolops.AssignStereochemistry(molecule)
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
        raise ValueError("File format not recognized. RDKit does not properly support mol2 files.")
    


    return [molecule]
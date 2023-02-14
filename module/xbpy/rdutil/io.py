
import logging

def read_molecules(path, recursive = True, *args, **kwargs):
    """Read molecules as RDK molecules from a single pdb, mol or sdf file, or from a directory /multiple directories of such files.
    
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
            print(p)
            for f in glob.glob(p):
                if os.path.isdir(f):
                    # If recursive we do a random walk through the directory tree
                    if (f not in seen) and recursive:
                        remaining_paths.add(f)
                else:
                    path_molecules.extend(_read_molecules_file(f, *args, **kwargs))
            seen.add(p)

        if len(path_molecules) == 0:
            logging.warning("No molecules found at path(s) {}. Make sure the file(s) exists.".format(path))
        else:
            molecules.extend(path_molecules)
    return molecules


def _read_molecules_file(path, *args, **kwargs):
    """Read molecules as RDK molecules from a single pdb, mol or sdf file.

    Args:
        path (str): Path to a file or directory. May be a list of paths. A path may also describe a glob pattern.

    Returns:
        list(RDKit.Mol): List of RDKit molecules.

    """

    import os
    import rdkit.Chem as Chem

    molecules = []
    if os.path.splitext(path)[1] == ".pdb":
        molecules = Chem.rdmolfiles.MolFromPDBFile(path, *args, **kwargs)
    elif os.path.splitext(path)[1] == ".mol":
        molecules = Chem.rdmolfiles.MolFromMolFile(path, *args, **kwargs)
    elif os.path.splitext(path)[1] == ".sdf":
        molecules = Chem.rdmolfiles.SDMolSupplier(path, *args, **kwargs)
    else:
        raise ValueError("File format not recognized. RDKit does not properly support mol2 files.")

    return molecules
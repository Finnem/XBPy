from scipy.spatial import cKDTree
from rdkit.Chem import AllChem as Chem
from ..rdutil import position
import numpy as np

class Receptor():
    """ Abstract class to represent an interacting receptor molecule.

    Parameters
    ----------
    mol : RDKit Mol
        The receptor molecule.

    Attributes
    ----------
    mol : RDKit Mol
        The receptor molecule.
    information : dict
        Information about the query molecule.
    query_positions : np.ndarray
        The positions of the atoms in the query molecule.
    kdtree : cKDTree
        A KDTree of the query positions.    

    """
    def __init__(self, mol):
        self.mol = mol
        self.information = []
        self.query_positions = []

        self._compute_query_positions_and_information(mol)
        self._construct_kd_tree()
        
    def _compute_query_positions_and_information(mol):
        """ Compute the positions of the atoms in the query molecule. """
        raise NotImplementedError("This method must be implemented in a subclass.")

    def _construct_kd_tree(self):
        """ Construct a KDTree from the query positions. """
        self.kdtree = cKDTree(self.query_positions)


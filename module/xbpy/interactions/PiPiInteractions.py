from .Base import Receptor
from scipy.spatial import cKDTree
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdmolops
from ..rdutil import is_planar
from ..rdutil import position
import numpy as np

class AromaticProximityReceptor(Receptor):
    """ Class to represent a receptor molecule that contains aromatic rings.
    Allows to detect aromatic-aromatic interactions between the receptor and a ligand via proximity of ring centers.

    Parameters
    ----------
    mol : RDKit Mol
        The receptor molecule.
    distance_threshold : float
        The distance threshold for pi pi proximity interactions.

    Attributes
    ----------  
    distance_threshold : float
        The distance threshold for pi pi proximity interactions.
    
    """

    def __init__(self, mol, distance_threshold = 5):
        self.distance_threshold = distance_threshold
        self.mol = mol
        super().__init__(mol)

    def _compute_query_positions_and_information(self, mol):
        """ Compute the positions of the atoms in the query molecule. 
        
        Parameters
        ----------
        mol : RDKit Mol
            The receptor molecule.

        """
        query_positions = []

        # compute rings
        try:
            Chem.SanitizeMol(mol)
        except:
            rdmolops.FastFindRings(mol)
        ri = mol.GetRingInfo()
        ring_indices = []
        for atom_indices in ri.AtomRings():
            atoms = [mol.GetAtomWithIdx(i) for i in atom_indices]
            if is_planar(atoms):
                query_positions.append(np.mean(position(atoms), axis = 0))
                ring_indices.append(atom_indices)

        self.query_positions = np.array(query_positions)
        self.information = {"ring_indices": np.array(ring_indices, dtype=object)}

    def detect_interactions(self, ligand):
        """ Detect interactions between the receptor and a ligand molecule.
        
        Parameters
        ----------
        ligand : RDKit Mol
            The ligand molecule.

        Returns
        -------
        list
            A list of interactions.

        """
        interactions = []
        ligand_positions = ligand.query_positions

        neighbors = self.kdtree.query_ball_point(ligand_positions, self.distance_threshold)
        other_indices = np.concatenate([np.full(len(indices), i, dtype=int) for i, indices in enumerate(neighbors)]).astype(int)
        this_indices = np.concatenate(neighbors).astype(int)

        if (len(this_indices) == 0) or (len(other_indices) == 0):
            return []
        
        other_atom_indices = ligand.information["ring_indices"][other_indices]
        this_atom_indices = self.information["ring_indices"][this_indices]

        other_positions = ligand.query_positions[other_indices]
        this_positions = self.query_positions[this_indices]

        return [(this_indices[i], other_indices[i], this_atom_indices[i], other_atom_indices[i], this_positions[i], other_positions[i]) for i in range(len(this_indices))]

class AromaticProximityLigand():
    """ Class to represent a ligand molecule that contains aromatic rings.
    Allows to detect aromatic-aromatic interactions between the receptor and a ligand via proximity of ring centers.

    Parameters
    ----------
    mol : RDKit Mol
        The receptor molecule.
    
    """

    def __init__(self, mol):
        self.mol = mol
        self._compute_query_positions_and_information(mol)

    def _compute_query_positions_and_information(self, mol):
        """ Compute the positions of the atoms in the query molecule. 
        
        Parameters
        ----------
        mol : RDKit Mol
            The receptor molecule.

        """
        query_positions = []

        # compute rings
        try:
            Chem.SanitizeMol(mol)
        except:
            rdmolops.FastFindRings(mol)
        ri = mol.GetRingInfo()
        ring_indices = []
        for atom_indices in ri.AtomRings():
            atoms = [mol.GetAtomWithIdx(i) for i in atom_indices]
            if is_planar(atoms):
                query_positions.append(np.mean(position(atoms), axis = 0))
                ring_indices.append(atom_indices)

        self.query_positions = np.array(query_positions)
        self.information = {"ring_indices": np.array(ring_indices, dtype=object)}



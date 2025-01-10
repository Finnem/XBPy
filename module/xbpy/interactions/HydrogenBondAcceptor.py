from .Base import Receptor
from scipy.spatial import cKDTree
from rdkit.Chem import AllChem as Chem
from ..rdutil import position
import numpy as np

class HydrogenBondAcceptorReceptor(Receptor):
    """ Class to represent a receptor molecule that can act as a hydrogen bond acceptor.

    Parameters
    ----------
    mol : RDKit Mol
        The receptor molecule.
    distance_threshold : float
        The distance threshold for hydrogen bond acceptor interactions.
    angle_threshold : float 
        The angle threshold for hydrogen bond acceptor interactions in degree.
    restricted_angle_threshold : float
        The angle threshold for hydrogen bond acceptor interactions in degree when the donor is restricted.

    Attributes
    ----------  
    query_mols : list
        A list of query molecules that represent hydrogen bond acceptors.
    distance_threshold : float
        The distance threshold for hydrogen bond acceptor interactions.
    angle_threshold : float
        The angle threshold for hydrogen bond acceptor interactions in degree.
    restricted_angle_threshold : float
        The angle threshold for hydrogen bond acceptor interactions in degree when the donor is restricted.
    """

    def __init__(self, mol, distance_threshold = 2.8, angle_threshold = 30, restricted_angle_threshold = 45):
        self.query_mols = [
            Chem.MolFromSmarts("[OX2;H1:1]"), # hydroxyl
            Chem.MolFromSmarts("[OX1;H0:1]"), # carbonyl
            Chem.MolFromSmarts("[OX2;H0:1]"), # ether
            Chem.MolFromSmarts("[NX3;H2:1]"), # amide
            Chem.MolFromSmarts("[NX3;H1:1]"), # amine
            Chem.MolFromSmarts("[NX2;H0:1]"), # imine
            Chem.MolFromSmarts("[NX3;H0:1]"), # amine
            Chem.MolFromSmarts("[SX2;H0:1]"), # thioether
            Chem.MolFromSmarts("[SX1;H1:1]"), # thiol
        ]
        self.distance_threshold = distance_threshold
        self.angle_threshold = angle_threshold
        self.restricted_angle_threshold = restricted_angle_threshold
        super().__init__(mol)

    def _compute_query_positions_and_information(self, mol):
        """ Compute the positions of the atoms in the query molecule. 
        
        Parameters
        ----------
        mol : RDKit Mol
            The receptor molecule.

        """
        query_positions = []
        restricted_directions = [[], [], []]
        restricted = [[], []]
        for query_mol in self.query_mols:
            for match in mol.GetSubstructMatches(query_mol):
                matched_atom = mol.GetAtomWithIdx(match[0])
                matched_position = position(matched_atom)
                neighbors = matched_atom.GetNeighbors()
                for i in range(1, len(neighbors) - 2, -1):
                    restricted[i].append(False)
                for i, neighbor in enumerate(matched_atom.GetNeighbors()):
                    neighbor_position = position(neighbor)
                    direction = neighbor_position - matched_position; direction /= np.linalg.norm(direction)
                    restricted_directions[i].append(direction)
                    if i > 0:
                        restricted[i - 1].append(True)
                query_positions.append(matched_position)
        restricted = [np.ones(len(query_positions), dtype=bool)] + restricted
        self.information = {
            "restricted_directions" : [np.array(l) for l in restricted_directions], 
            "restricted" : np.array(restricted)
        }
        self.query_positions = np.array(query_positions)

    def detect_interactions(self, ligand):
        """ Detects hydrogen bond acceptor interactions between the receptor and a ligand. 
        
        Parameters
        ----------
        ligand : HydrogenBondDonorLigand
            The ligand molecule.

        Returns
        -------
        list
            A list of hydrogen bond acceptor interactions.
        """

        # detect if hydrogens should be ignored
        ignore_Hs = ligand.ignore_Hs
        if ignore_Hs:
            distance_threshold = self.distance_threshold + 1.1
        else:
            distance_threshold = self.distance_threshold

        # Find all systems within the distance threshold
        acceptor_indices = self.kdtree.query_ball_point(ligand.query_positions, self.distance_threshold)
        donor_indices = np.concatenate([np.full(len(indices), i, dtype=int) for i, indices in enumerate(acceptor_indices)])
        acceptor_indices = np.concatenate(acceptor_indices).astype(int)

        # compute the angle between the hydrogen bond acceptor and the hydrogen bond donor
        differences = self.query_positions[acceptor_indices] - ligand.query_positions[donor_indices]
        distances = np.linalg.norm(differences, axis = 1)
        acceptor_directions = differences / distances[:, np.newaxis]
        acceptor_positions = self.query_positions[acceptor_indices]
        donor_positions = ligand.query_positions[donor_indices]

        if ignore_Hs:
            donor_directions = acceptor_directions
        else:
            donor_directions = ligand.information["directions"][donor_indices]
        angles = np.arccos(np.clip(np.sum(acceptor_directions * donor_directions, axis = 1), -1.0, 1.0)) * 180 / np.pi

        # filter out interactions that do not meet the angle threshold
        mask = angles < self.angle_threshold
        acceptor_indices = acceptor_indices[mask]
        acceptor_directions = acceptor_directions[mask]
        donor_indices = donor_indices[mask]
        angles = angles[mask]
        distances = distances[mask]
        acceptor_positions = acceptor_positions[mask]
        donor_positions = donor_positions[mask]

        # check first neighbor restriction
        
        for i in range(3):
            # determine restricted direction mapping
            restricted_mask = self.information["restricted"][i]
            restricted_directions = self.information["restricted_directions"][i]
            inverse_restricted_directions = np.zeros((len(restricted_mask), 3))
            inverse_restricted_directions[restricted_mask] = restricted_directions

            # determined used restrictions
            used_restricted = restricted_mask[acceptor_indices]
            used_restricted_directions = inverse_restricted_directions[acceptor_indices]
            used_restricted_directions = used_restricted_directions[used_restricted]

            # compute the angle
            restricted_angles = np.arccos(np.clip(np.sum(-acceptor_directions[used_restricted] * used_restricted_directions, axis = 1), -1.0, 1.0)) * 180 / np.pi
            used_restricted[used_restricted] = restricted_angles < self.restricted_angle_threshold
            mask = ~used_restricted

            acceptor_indices = acceptor_indices[mask]
            acceptor_directions = acceptor_directions[mask]
            donor_indices = donor_indices[mask]
            angles = angles[mask]
            distances = distances[mask]
            acceptor_positions = acceptor_positions[mask]
            donor_positions = donor_positions[mask]

        return list(zip(acceptor_indices, donor_indices, acceptor_positions, donor_positions, angles, distances))

class HydrogenBondDonorLigand():
    """ Class to represent a ligand molecule that can act as a hydrogen bond donor.

    Parameters
    ----------
    mol : RDKit Mol
        The ligand molecule.
    
    Attributes
    ----------  
    query_mols : list
        A list of query molecules that represent hydrogen bond donors.  
    mol : RDKit Mol
        The ligand molecule.
    information : dict
        Information about the query molecule.
    query_positions : np.ndarray
        The positions of the atoms in the query molecule.
    ignore_Hs : bool
        Whether to ignore hydrogens in the query molecule. This means that ideal hydrogen bond donors are assumed.

    """

    def __init__(self, mol, ignore_Hs = False):
        self.query_mols = [
            Chem.MolFromSmarts("[OX2;H1][H:1]"), # hydroxyl
            Chem.MolFromSmarts("[NX2;H1][H:1]"), # amine
            Chem.MolFromSmarts("[NX3;H1][H:1]"), # amine
            Chem.MolFromSmarts("[NX3;H2][H:1]"), # amide
            Chem.MolFromSmarts("[SX2;H1][H:1]"), # thiol
        ]
        self.mol = mol
        self.ignore_Hs = ignore_Hs
        self.information = []
        self.query_positions = []

        self._compute_query_positions_and_information(mol)
        
    def _compute_query_positions_and_information(self, mol):
        """ Compute the positions of the atoms in the query molecule. """
        query_positions = []
        directions = []
        for query_mol in self.query_mols:
            for match in mol.GetSubstructMatches(query_mol):
                if self.ignore_Hs:
                    matched_atom = mol.GetAtomWithIdx(match[0])
                else:
                    matched_atom = mol.GetAtomWithIdx(match[1]) # not sure why first group is still captured
                matched_position = position(matched_atom)
                neighbors = matched_atom.GetNeighbors()
                neighbor_position = position(neighbors[0])
                direction = matched_position - neighbor_position; direction /= np.linalg.norm(direction)
                directions.append(direction)
                query_positions.append(matched_position)
        self.information = {
            "directions" : np.array(directions)
        }
        self.query_positions = np.array(query_positions)

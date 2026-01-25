from xbpy import rdutil
from xbpy.orcautil import fragment_molecule
import rdkit.Chem as Chem
import numpy as np
import pymolviz as pmv
import logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(asctime)s - %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def test_fragmentation():
    #mol = next(rdutil.read_molecules("larger_test_case.xyz"))
    #rdutil.write_molecules([mol], "larger_test_case_connected.sdf")
    mol = next(rdutil.read_molecules("larger_test_case_connected.sdf"))
    frags, to_replace = fragment_molecule(mol, 9329, verbose = True)
    all_indices = set()
    for frag in frags:
        all_indices.update(frag)
    remaining_mol = Chem.RWMol(rdutil.keep_atoms(mol, all_indices))
    index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted(list(all_indices)))}
    # replace to_replace atoms with hydrogens at distance 1.1 and bond order 1
    for atom_pair in to_replace:
        atom1 = mol.GetAtomWithIdx(atom_pair[0])
        atom2 = mol.GetAtomWithIdx(atom_pair[1])
        remaining_mol.AddAtom(Chem.Atom(1))
        new_idx = remaining_mol.GetNumAtoms() - 1
        remaining_mol.AddBond(index_map[atom1.GetIdx()], new_idx, Chem.rdchem.BondType.SINGLE)
        pos1 = rdutil.position(atom1)
        pos2 = rdutil.position(atom2)
        direction = pos2 - pos1
        direction /= np.linalg.norm(direction)
        new_position = pos1 + direction * 1.1
        remaining_mol.GetConformer().SetAtomPosition(new_idx, new_position)
    rdutil.write_molecules([remaining_mol], "larger_test_case_remaining.sdf")
    i = 0
    for frag in frags:
        pmv.Points(rdutil.position(mol)[frag], name = f"frag_{i}").write(f"frag_{i}.py")
        i += 1

if __name__ == "__main__":
    test_fragmentation()
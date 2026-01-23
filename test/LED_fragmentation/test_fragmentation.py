from xbpy import rdutil
from xbpy.orcautil import fragment_molecule
import pymolviz as pmv
import logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(asctime)s - %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def test_fragmentation():
    mol = next(rdutil.read_molecules("larger_test_case.xyz", verbose = True))
    frags = fragment_molecule(mol, 1519)
    i = 0
    for frag in frags:
        pmv.Points(rdutil.position(mol)[frag]).write(f"test/LED_fragmentation/frag_{i}.py")
        i += 1

if __name__ == "__main__":
    test_fragmentation()

from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


ligand_aromaticity_centers = [
        
COLOR,0.21044753832183283,0.6773105080456748,0.6433941168468681,1.0,SPHERE,11.516666666666666,11.029499999999999,7.591333333333334,0.3

            ]
cmd.load_cgo(ligand_aromaticity_centers, "ligand_aromaticity_centers", state=1)
        

cmd.set("cgo_transparency", 0, "ligand_aromaticity_centers")
    

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

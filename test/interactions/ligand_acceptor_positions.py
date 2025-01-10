
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


ligand_acceptor_positions = [
        
COLOR,1.0,0.0,0.0,1.0,SPHERE,14.002,20.149,8.594,0.3,COLOR,1.0,0.0,0.0,1.0,SPHERE,10.261,15.844,13.217,0.3

            ]
cmd.load_cgo(ligand_acceptor_positions, "ligand_acceptor_positions", state=1)
        

cmd.set("cgo_transparency", 0, "ligand_acceptor_positions")
    

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

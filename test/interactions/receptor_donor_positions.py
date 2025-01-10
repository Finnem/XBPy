
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


receptor_donor_positions = [
        
COLOR,0.0,0.0,1.0,1.0,SPHERE,15.832,21.137,9.631,0.3,COLOR,0.0,0.0,1.0,1.0,SPHERE,11.791,14.12,14.707,0.3

            ]
cmd.load_cgo(receptor_donor_positions, "receptor_donor_positions", state=1)
        

cmd.set("cgo_transparency", 0, "receptor_donor_positions")
    

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

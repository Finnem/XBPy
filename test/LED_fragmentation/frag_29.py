
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_29 = [
        
COLOR,0.86,0.33999999999999997,0.46479999999999977,1.0,SPHERE,45.849998,73.419998,53.080002,0.3,COLOR,0.86,0.33999999999999997,0.46479999999999977,1.0,SPHERE,45.220001,73.980003,52.619999,0.3,COLOR,0.86,0.33999999999999997,0.46479999999999977,1.0,SPHERE,46.709999,73.779999,52.849998,0.3

            ]
cmd.load_cgo(frag_29, "frag_29", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

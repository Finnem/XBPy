
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_3 = [
        
COLOR,0.41279999999999994,0.86,0.33999999999999997,1.0,SPHERE,39.250477,66.438255,53.605328,0.3,COLOR,0.41279999999999994,0.86,0.33999999999999997,1.0,SPHERE,40.15057,66.966415,53.946297,0.3,COLOR,0.41279999999999994,0.86,0.33999999999999997,1.0,SPHERE,39.531754,65.405632,53.353367,0.3,COLOR,0.41279999999999994,0.86,0.33999999999999997,1.0,SPHERE,38.884205,66.934364,52.694275,0.3

            ]
cmd.load_cgo(frag_3, "frag_3", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

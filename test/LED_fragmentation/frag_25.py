
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_25 = [
        
COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,41.330002,63.310001,47.830002,0.3,COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,40.490002,63.400002,47.369999,0.3,COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,41.119999,62.790001,48.610001,0.3

            ]
cmd.load_cgo(frag_25, "frag_25", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

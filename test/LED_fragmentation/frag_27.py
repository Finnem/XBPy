
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_27 = [
        
COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,40.34,70.849998,45.080002,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,40.860001,70.639999,44.310001,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,40.880001,71.459999,45.580002,0.3

            ]
cmd.load_cgo(frag_27, "frag_27", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

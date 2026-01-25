
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_34 = [
        
COLOR,0.33999999999999997,0.86,0.8132000000000001,1.0,SPHERE,45.669998,71.540001,57.0,0.3,COLOR,0.33999999999999997,0.86,0.8132000000000001,1.0,SPHERE,46.299999,71.480003,57.720001,0.3,COLOR,0.33999999999999997,0.86,0.8132000000000001,1.0,SPHERE,45.439999,70.629997,56.82,0.3

            ]
cmd.load_cgo(frag_34, "frag_34", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

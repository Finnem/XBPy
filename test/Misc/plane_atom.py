
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


Points_0 = [
        
COLOR,1.0,0.0,0.0,1.0,SPHERE,1.192,-0.465,0.67,0.3

            ]
cmd.load_cgo(Points_0, "Points_0", state=1)
cmd.set("cgo_transparency", 0, "Points_0")
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

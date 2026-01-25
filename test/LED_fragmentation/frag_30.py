
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_30 = [
        
COLOR,0.86,0.6052,0.33999999999999997,1.0,SPHERE,42.470001,72.089996,46.639999,0.3,COLOR,0.86,0.6052,0.33999999999999997,1.0,SPHERE,42.689999,72.449997,47.5,0.3,COLOR,0.86,0.6052,0.33999999999999997,1.0,SPHERE,42.779999,71.18,46.68,0.3

            ]
cmd.load_cgo(frag_30, "frag_30", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

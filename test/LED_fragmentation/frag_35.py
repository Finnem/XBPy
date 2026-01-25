
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_35 = [
        
COLOR,0.33999999999999997,0.5947999999999997,0.86,1.0,SPHERE,44.27,73.959999,55.279999,0.3,COLOR,0.33999999999999997,0.5947999999999997,0.86,1.0,SPHERE,43.509998,73.370003,55.259998,0.3,COLOR,0.33999999999999997,0.5947999999999997,0.86,1.0,SPHERE,44.810001,73.669998,54.549999,0.3

            ]
cmd.load_cgo(frag_35, "frag_35", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

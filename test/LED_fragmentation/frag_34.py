
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_34 = [
        
COLOR,0.33999999999999997,0.86,0.8132000000000001,1.0,SPHERE,55.43,97.449997,93.209999,0.3,COLOR,0.33999999999999997,0.86,0.8132000000000001,1.0,SPHERE,54.93,98.220001,93.489998,0.3,COLOR,0.33999999999999997,0.86,0.8132000000000001,1.0,SPHERE,56.139999,97.389999,93.849998,0.3

            ]
cmd.load_cgo(frag_34, "frag_34", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

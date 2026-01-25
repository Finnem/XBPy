
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_5 = [
        
COLOR,0.33999999999999997,0.8287999999999999,0.86,1.0,SPHERE,38.217773,64.862656,58.175972,0.3,COLOR,0.33999999999999997,0.8287999999999999,0.86,1.0,SPHERE,39.165798,64.335068,58.041405,0.3,COLOR,0.33999999999999997,0.8287999999999999,0.86,1.0,SPHERE,37.39711901585309,64.13909520601936,58.062012289058316,0.3,COLOR,0.33999999999999997,0.8287999999999999,0.86,1.0,SPHERE,38.17530562338896,65.31571849922437,59.17743637466359,0.3

            ]
cmd.load_cgo(frag_5, "frag_5", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

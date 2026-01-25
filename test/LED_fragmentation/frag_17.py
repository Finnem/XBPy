
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_17 = [
        
COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,49.76701,69.963814,57.336918,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,50.653767,69.396431,57.064297,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,50.068878,70.689468,58.093807,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,49.03961016893268,69.27113484760456,57.785343108128394,0.3

            ]
cmd.load_cgo(frag_17, "frag_17", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

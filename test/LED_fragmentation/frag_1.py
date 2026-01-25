
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_1 = [
        
COLOR,0.86,0.6832,0.33999999999999997,1.0,SPHERE,40.314556,74.540749,48.142361,0.3,COLOR,0.86,0.6832,0.33999999999999997,1.0,SPHERE,40.388584,75.078316,47.186405,0.3,COLOR,0.86,0.6832,0.33999999999999997,1.0,SPHERE,41.293736,74.543777,48.634186,0.3,COLOR,0.86,0.6832,0.33999999999999997,1.0,SPHERE,40.011597,73.501129,47.947044,0.3

            ]
cmd.load_cgo(frag_1, "frag_1", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

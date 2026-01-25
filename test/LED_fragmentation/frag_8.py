
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_8 = [
        
COLOR,0.7871999999999999,0.33999999999999997,0.86,1.0,SPHERE,46.563175,66.089203,54.819458,0.3,COLOR,0.7871999999999999,0.33999999999999997,0.86,1.0,SPHERE,47.629128,66.331978,54.827255,0.3,COLOR,0.7871999999999999,0.33999999999999997,0.86,1.0,SPHERE,46.071278,66.830589,54.175533,0.3,COLOR,0.7871999999999999,0.33999999999999997,0.86,1.0,SPHERE,46.410568,64.422035,54.079502,0.3,COLOR,0.7871999999999999,0.33999999999999997,0.86,1.0,SPHERE,45.061634,64.371216,54.118568,0.3

            ]
cmd.load_cgo(frag_8, "frag_8", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

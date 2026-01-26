
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_26 = [
        
COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,43.139999,72.940002,49.830002,0.3,COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,42.75,73.620003,50.389999,0.3,COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,44.040001,72.860001,50.130001,0.3

            ]
cmd.load_cgo(frag_26, "frag_26", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

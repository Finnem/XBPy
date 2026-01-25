
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_16 = [
        
COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,49.28257,70.747047,56.09977,0.3,COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,49.915619,71.746811,55.773685,0.3,COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,48.156452,70.334068,55.440475,0.3,COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,47.719261,69.472336,55.7234,0.3,COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,47.702625,70.979523,54.251835,0.3,COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,47.824589,72.074539,54.304871,0.3,COLOR,0.33999999999999997,0.36079999999999973,0.86,1.0,SPHERE,46.614952,70.828552,54.137451,0.3

            ]
cmd.load_cgo(frag_16, "frag_16", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

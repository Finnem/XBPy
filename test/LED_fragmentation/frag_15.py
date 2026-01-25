
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_15 = [
        
COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,43.834141,65.20858,56.985504,0.3,COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,44.373318,64.50528,57.836071,0.3,COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,44.518562,65.957222,56.110271,0.3,COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,44.015736,66.416573,55.32999,0.3,COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,45.959911,66.090469,56.222122,0.3,COLOR,0.33999999999999997,0.6727999999999997,0.86,1.0,SPHERE,46.3288,65.228729,56.794334,0.3

            ]
cmd.load_cgo(frag_15, "frag_15", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

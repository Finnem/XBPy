
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_9 = [
        
COLOR,0.86,0.33999999999999997,0.6207999999999999,1.0,SPHERE,44.27243,64.733032,47.646645,0.3,COLOR,0.86,0.33999999999999997,0.6207999999999999,1.0,SPHERE,44.737289,64.401886,48.585426,0.3,COLOR,0.86,0.33999999999999997,0.6207999999999999,1.0,SPHERE,43.193962,64.83902,47.821648,0.3,COLOR,0.86,0.33999999999999997,0.6207999999999999,1.0,SPHERE,44.862144,66.422241,47.254711,0.3,COLOR,0.86,0.33999999999999997,0.6207999999999999,1.0,SPHERE,46.134197,66.044235,46.958115,0.3

            ]
cmd.load_cgo(frag_9, "frag_9", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()


from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_14 = [
        
COLOR,0.33999999999999997,0.86,0.7352000000000001,1.0,SPHERE,34.939167,61.688984,53.301521,0.3,COLOR,0.33999999999999997,0.86,0.7352000000000001,1.0,SPHERE,35.06506,62.77002,53.386475,0.3,COLOR,0.33999999999999997,0.86,0.7352000000000001,1.0,SPHERE,34.23596069063097,61.3514143911314,54.077116026822566,0.3,COLOR,0.33999999999999997,0.86,0.7352000000000001,1.0,SPHERE,34.52599297123164,61.44877202640522,52.31077106702346,0.3

            ]
cmd.load_cgo(frag_14, "frag_14", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

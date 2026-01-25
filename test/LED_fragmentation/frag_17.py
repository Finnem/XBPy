
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_17 = [
        
COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,41.879005,63.172955,55.793739,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,42.189026,63.614571,54.682209,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,41.796532,63.945297,56.886738,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,41.706802,63.474625,57.793129,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,42.3158,65.290924,56.873405,0.3,COLOR,0.6311999999999998,0.33999999999999997,0.86,1.0,SPHERE,42.030914,65.758446,55.921997,0.3

            ]
cmd.load_cgo(frag_17, "frag_17", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

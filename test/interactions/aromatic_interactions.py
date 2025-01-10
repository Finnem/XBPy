
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


aromatic_interactions = [
        
LINEWIDTH,1,BEGIN,LINES,COLOR,0.8004936186423958,0.47703363533737203,0.9579547196007522,VERTEX,9.9535,9.889000000000001,10.734666666666667,COLOR,0.8004936186423958,0.47703363533737203,0.9579547196007522,VERTEX,11.516666666666666,11.029499999999999,7.591333333333334,END

            ]
cmd.load_cgo(aromatic_interactions, "aromatic_interactions", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

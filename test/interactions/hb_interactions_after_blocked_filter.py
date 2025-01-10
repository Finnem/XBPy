
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


after_blocked_filter = [
        
LINEWIDTH,1,BEGIN,LINES,COLOR,0.9677975592919913,0.44127456009157356,0.5358103155058701,VERTEX,11.791,14.12,14.707,COLOR,0.9677975592919913,0.44127456009157356,0.5358103155058701,VERTEX,10.261,15.844,13.217,END

            ]
cmd.load_cgo(after_blocked_filter, "after_blocked_filter", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

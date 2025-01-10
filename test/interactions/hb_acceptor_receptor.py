
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


Detected_Acceptors = [
        
COLOR,0.21044753832183283,0.6773105080456748,0.6433941168468681,1.0,SPHERE,14.03,15.623,4.918,0.3,COLOR,0.21044753832183283,0.6773105080456748,0.6433941168468681,1.0,SPHERE,32.177,6.916,4.524,0.3,COLOR,0.21044753832183283,0.6773105080456748,0.6433941168468681,1.0,SPHERE,16.472,17.78,6.465,0.3,COLOR,0.21044753832183283,0.6773105080456748,0.6433941168468681,1.0,SPHERE,27.754,6.87,3.017,0.3

            ]
cmd.load_cgo(Detected_Acceptors, "Detected_Acceptors", state=1)
        

cmd.set("cgo_transparency", 0, "Detected_Acceptors")
    

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

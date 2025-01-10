
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


Detected_Donors = [
        
COLOR,0.46810256823426105,0.6699492535792404,0.1928958739904499,1.0,SPHERE,13.612,14.265,4.498,0.3,COLOR,0.46810256823426105,0.6699492535792404,0.1928958739904499,1.0,SPHERE,13.612,14.265,4.498,0.3,COLOR,0.46810256823426105,0.6699492535792404,0.1928958739904499,1.0,SPHERE,16.142,19.978,5.947,0.3,COLOR,0.46810256823426105,0.6699492535792404,0.1928958739904499,1.0,SPHERE,16.142,19.978,5.947,0.3

            ]
cmd.load_cgo(Detected_Donors, "Detected_Donors", state=1)
        

cmd.set("cgo_transparency", 0, "Detected_Donors")
    

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

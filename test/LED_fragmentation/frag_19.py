
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_19 = [
        
COLOR,0.86,0.33999999999999997,0.46479999999999977,1.0,SPHERE,48.480064,71.507591,44.958218,0.3,COLOR,0.86,0.33999999999999997,0.46479999999999977,1.0,SPHERE,48.621559,72.024681,45.908077,0.3,COLOR,0.86,0.33999999999999997,0.46479999999999977,1.0,SPHERE,47.891764116968176,72.15062516822157,44.287091852301714,0.3,COLOR,0.86,0.33999999999999997,0.46479999999999977,1.0,SPHERE,49.46393244013662,71.30460670260184,44.5101105714598,0.3

            ]
cmd.load_cgo(frag_19, "frag_19", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

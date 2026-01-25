
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_20 = [
        
COLOR,0.86,0.5272,0.33999999999999997,1.0,SPHERE,40.180103,74.410065,51.259895,0.3,COLOR,0.86,0.5272,0.33999999999999997,1.0,SPHERE,41.346867,74.732018,51.010464,0.3,COLOR,0.86,0.5272,0.33999999999999997,1.0,SPHERE,39.175907,74.562065,50.371181,0.3,COLOR,0.86,0.5272,0.33999999999999997,1.0,SPHERE,38.23978,74.293694,50.66428,0.3,COLOR,0.86,0.5272,0.33999999999999997,1.0,SPHERE,39.29953,75.227631,49.082756,0.3,COLOR,0.86,0.5272,0.33999999999999997,1.0,SPHERE,39.622135,76.276001,49.229546,0.3

            ]
cmd.load_cgo(frag_20, "frag_20", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

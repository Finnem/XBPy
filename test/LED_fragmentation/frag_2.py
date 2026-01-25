
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


frag_2 = [
        
COLOR,0.7247999999999999,0.86,0.33999999999999997,1.0,SPHERE,36.815067,75.623322,54.242458,0.3,COLOR,0.7247999999999999,0.86,0.33999999999999997,1.0,SPHERE,36.153824,74.814049,53.937515,0.3,COLOR,0.7247999999999999,0.86,0.33999999999999997,1.0,SPHERE,36.56231850280958,76.51873695608857,53.65564435532148,0.3,COLOR,0.7247999999999999,0.86,0.33999999999999997,1.0,SPHERE,36.664041281639236,75.82469655610456,55.313270551545465,0.3

            ]
cmd.load_cgo(frag_2, "frag_2", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

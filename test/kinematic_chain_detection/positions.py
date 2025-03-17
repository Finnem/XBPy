
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


Positions = [
        
COLOR,0.46810256823426105,0.6699492535792404,0.1928958739904499,1.0,SPHERE,7.6101510697585635,3.897229052943425,0.5566446532772105,0.3,COLOR,0.46810256823426105,0.6699492535792404,0.1928958739904499,1.0,SPHERE,6.686121465155587,4.5116686563908655,1.0926290406735406,0.3

            ]
cmd.load_cgo(Positions, "Positions", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

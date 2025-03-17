
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


Arrows_5 = [
        
ALPHA,1,CONE,9.676747450603244,18.160409741684603,5.975261918555509,11.363461862650812,18.48282743542115,5.936515479638878,0.05,0.05,0.7632105624545802,0.5838460616396939,0.19465686802007026,0.7632105624545802,0.5838460616396939,0.19465686802007026,1.0,0.0,ALPHA,1,CONE,11.363461862650812,18.48282743542115,5.936515479638878,11.9257,18.5903,5.9236,0.08090000000000001,0.0,0.7632105624545802,0.5838460616396939,0.19465686802007026,0.7632105624545802,0.5838460616396939,0.19465686802007026,0.0,0.0

            ]
cmd.load_cgo(Arrows_5, "Arrows_5", state=1)
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

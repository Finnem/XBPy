
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


Arrows_0 = [
        
CYLINDER,0.0,0.0,0.0,-0.3755760112233604,-0.4904095538534062,-0.425371754213715,0.3,1.0,0.0,0.0,1.0,0.0,0.0,CYLINDER,0.0,0.0,0.0,0.0,-0.4904095538534062,-0.425371754213715,0.3,0.01568627450980431,0.4940868896578237,0.0,0.0,0.5019607843137255,0.0,CONE,-0.3755760112233604,-0.4904095538534062,-0.425371754213715,-0.5007680149644805,-0.6538794051378749,-0.5671623389516199,0.4854,0.0,1.0,0.0,0.0,1.0,0.0,0.0,1.0,1.0,CONE,0.0,-0.4904095538534062,-0.425371754213715,0.0,-0.6538794051378749,-0.5671623389516199,0.4854,0.0,0.0,0.5019607843137255,0.0,0.0,0.5019607843137255,0.0,1.0,1.0

            ]
cmd.load_cgo(Arrows_0, "Arrows_0", state=1)
cmd.set("cgo_transparency", 0, "Arrows_0")
        

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

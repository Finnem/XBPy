
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


problem_cluster_0 = [
        
COLOR,0.9677975592919913,0.44127456009157356,0.5358103155058701,1.0,SPHERE,-0.713,0.0,1.154,0.3,COLOR,0.9677975592919913,0.44127456009157356,0.5358103155058701,1.0,SPHERE,-0.039,0.05,2.435,0.3,COLOR,0.9677975592919913,0.44127456009157356,0.5358103155058701,1.0,SPHERE,0.402,1.03,2.616,0.3,COLOR,0.9677975592919913,0.44127456009157356,0.5358103155058701,1.0,SPHERE,0.255,1.883,1.882,0.3

            ]
cmd.load_cgo(problem_cluster_0, "problem_cluster_0", state=1)
        

cmd.set("cgo_transparency", 0, "problem_cluster_0")
    

problem_clusters = cmd.group("problem_clusters")
cmd.group("problem_clusters", "open")

cmd.group("problem_clusters", "problem_cluster_0", "add")

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

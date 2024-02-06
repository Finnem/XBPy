
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))

cmd.pseudoatom("mol1_history", label="1.0", pos = [-8.785, 3.378, 0.725], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [0.16, 1.309, 0.392], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-9.493, 4.534, 0.726], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-7.558, 3.378, 0.726], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-9.591, 2.103, 0.732], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-8.812, 5.812, 0.687], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-1.975, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-1.292, 1.309, 0.298], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-1.292, 1.309, -2.13], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [0.101, 1.309, -2.121], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [0.8, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-3.745, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-10.493, 4.495, 0.675], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-9.534, 6.601, 0.871], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-8.043, 5.839, 1.454], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-8.333, 5.976, -0.276], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-10.664, 2.274, 0.762], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-9.295, 1.514, 1.596], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-9.342, 1.532, -0.159], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-1.841, 1.309, 1.228], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-1.841, 1.309, -3.06], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [0.636, 1.309, -3.061], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [1.881, 1.309, -0.917], state = 1)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-8.785, 3.378, 0.725], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="2.0", pos = [0.16, 1.309, 0.392], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-9.493, 4.534, 0.726], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-7.558, 3.378, 0.726], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="4.0", pos = [-9.591, 2.103, 0.732], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="4.0", pos = [-8.812, 5.812, 0.687], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-1.975, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-1.292, 1.309, 0.298], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-1.292, 1.309, -2.13], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [0.101, 1.309, -2.121], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [0.8, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-3.745, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-10.493, 4.495, 0.675], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-9.534, 6.601, 0.871], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-8.043, 5.839, 1.454], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-8.333, 5.976, -0.276], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-10.664, 2.274, 0.762], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-9.295, 1.514, 1.596], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-9.342, 1.532, -0.159], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-1.841, 1.309, 1.228], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [-1.841, 1.309, -3.06], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [0.636, 1.309, -3.061], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="1.0", pos = [1.881, 1.309, -0.917], state = 2)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="8.0", pos = [-8.785, 3.378, 0.725], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [0.16, 1.309, 0.392], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="8.0", pos = [-9.493, 4.534, 0.726], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-7.558, 3.378, 0.726], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-9.591, 2.103, 0.732], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-8.812, 5.812, 0.687], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="7.0", pos = [-1.975, 1.309, -0.916], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-1.292, 1.309, 0.298], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="7.0", pos = [-1.292, 1.309, -2.13], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="7.0", pos = [0.101, 1.309, -2.121], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [0.8, 1.309, -0.916], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-3.745, 1.309, -0.916], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-10.493, 4.495, 0.675], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="4.0", pos = [-9.534, 6.601, 0.871], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="4.0", pos = [-8.043, 5.839, 1.454], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="4.0", pos = [-8.333, 5.976, -0.276], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="4.0", pos = [-10.664, 2.274, 0.762], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="4.0", pos = [-9.295, 1.514, 1.596], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="4.0", pos = [-9.342, 1.532, -0.159], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-1.841, 1.309, 1.228], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [-1.841, 1.309, -3.06], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [0.636, 1.309, -3.061], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="3.0", pos = [1.881, 1.309, -0.917], state = 3)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="17.0", pos = [-8.785, 3.378, 0.725], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="12.0", pos = [0.16, 1.309, 0.392], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="17.0", pos = [-9.493, 4.534, 0.726], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="8.0", pos = [-7.558, 3.378, 0.726], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="20.0", pos = [-9.591, 2.103, 0.732], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="20.0", pos = [-8.812, 5.812, 0.687], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="16.0", pos = [-1.975, 1.309, -0.916], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="16.0", pos = [-1.292, 1.309, 0.298], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="17.0", pos = [-1.292, 1.309, -2.13], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="16.0", pos = [0.101, 1.309, -2.121], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="16.0", pos = [0.8, 1.309, -0.916], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="7.0", pos = [-3.745, 1.309, -0.916], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="8.0", pos = [-10.493, 4.495, 0.675], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-9.534, 6.601, 0.871], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-8.043, 5.839, 1.454], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-8.333, 5.976, -0.276], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-10.664, 2.274, 0.762], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-9.295, 1.514, 1.596], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-9.342, 1.532, -0.159], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [-1.841, 1.309, 1.228], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="7.0", pos = [-1.841, 1.309, -3.06], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="7.0", pos = [0.636, 1.309, -3.061], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="6.0", pos = [1.881, 1.309, -0.917], state = 4)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="45.0", pos = [-8.785, 3.378, 0.725], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="32.0", pos = [0.16, 1.309, 0.392], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="45.0", pos = [-9.493, 4.534, 0.726], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="17.0", pos = [-7.558, 3.378, 0.726], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="35.0", pos = [-9.591, 2.103, 0.732], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="35.0", pos = [-8.812, 5.812, 0.687], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="40.0", pos = [-1.975, 1.309, -0.916], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="34.0", pos = [-1.292, 1.309, 0.298], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="39.0", pos = [-1.292, 1.309, -2.13], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="40.0", pos = [0.101, 1.309, -2.121], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="34.0", pos = [0.8, 1.309, -0.916], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="16.0", pos = [-3.745, 1.309, -0.916], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="17.0", pos = [-10.493, 4.495, 0.675], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="20.0", pos = [-9.534, 6.601, 0.871], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="20.0", pos = [-8.043, 5.839, 1.454], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="20.0", pos = [-8.333, 5.976, -0.276], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="20.0", pos = [-10.664, 2.274, 0.762], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="20.0", pos = [-9.295, 1.514, 1.596], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="20.0", pos = [-9.342, 1.532, -0.159], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="16.0", pos = [-1.841, 1.309, 1.228], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="17.0", pos = [-1.841, 1.309, -3.06], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="16.0", pos = [0.636, 1.309, -3.061], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="16.0", pos = [1.881, 1.309, -0.917], state = 5)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="97.0", pos = [-8.785, 3.378, 0.725], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="68.0", pos = [0.16, 1.309, 0.392], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="97.0", pos = [-9.493, 4.534, 0.726], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="45.0", pos = [-7.558, 3.378, 0.726], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="105.0", pos = [-9.591, 2.103, 0.732], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="105.0", pos = [-8.812, 5.812, 0.687], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="89.0", pos = [-1.975, 1.309, -0.916], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="88.0", pos = [-1.292, 1.309, 0.298], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="97.0", pos = [-1.292, 1.309, -2.13], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="89.0", pos = [0.101, 1.309, -2.121], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="88.0", pos = [0.8, 1.309, -0.916], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="40.0", pos = [-3.745, 1.309, -0.916], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="45.0", pos = [-10.493, 4.495, 0.675], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="35.0", pos = [-9.534, 6.601, 0.871], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="35.0", pos = [-8.043, 5.839, 1.454], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="35.0", pos = [-8.333, 5.976, -0.276], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="35.0", pos = [-10.664, 2.274, 0.762], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="35.0", pos = [-9.295, 1.514, 1.596], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="35.0", pos = [-9.342, 1.532, -0.159], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="34.0", pos = [-1.841, 1.309, 1.228], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="39.0", pos = [-1.841, 1.309, -3.06], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="40.0", pos = [0.636, 1.309, -3.061], state = 6)

cmd.set("label_size", 24, "mol1_history")

cmd.pseudoatom("mol1_history", label="34.0", pos = [1.881, 1.309, -0.917], state = 6)

cmd.set("label_size", 24, "mol1_history")


history1_group = cmd.group("history1_group")
cmd.group("history1_group", "open")

cmd.group("history1_group", "mol1_history", "add"),
cmd.group("history1_group", "mol1_history", "add"),
cmd.group("history1_group", "mol1_history", "add"),
cmd.group("history1_group", "mol1_history", "add"),
cmd.group("history1_group", "mol1_history", "add"),
cmd.group("history1_group", "mol1_history", "add")

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

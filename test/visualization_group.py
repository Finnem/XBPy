
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-7.558, 3.378, 0.726], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-3.745, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-10.493, 4.495, 0.675], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-9.534, 6.601, 0.871], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-8.043, 5.839, 1.454], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-8.333, 5.976, -0.276], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-10.664, 2.274, 0.762], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-9.295, 1.514, 1.596], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-9.342, 1.532, -0.159], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-1.841, 1.309, 1.228], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [-1.841, 1.309, -3.06], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [0.636, 1.309, -3.061], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1))", pos = [1.881, 1.309, -0.917], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 2))", pos = [0.16, 1.309, 0.392], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3))", pos = [-8.785, 3.378, 0.725], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3))", pos = [-9.493, 4.534, 0.726], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3))", pos = [-1.975, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3))", pos = [-1.292, 1.309, 0.298], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3))", pos = [-1.292, 1.309, -2.13], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3))", pos = [0.101, 1.309, -2.121], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3))", pos = [0.8, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 4))", pos = [-9.591, 2.103, 0.732], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 4))", pos = [-8.812, 5.812, 0.687], state = 1)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-6.93, 7.624, -1.575], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-1.106, 6.584, -3.75], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-1.106, 6.584, 0.538], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-5.779, 4.661, -2.409], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-5.732, 4.643, -0.654], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-7.101, 5.403, -1.488], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [1.371, 6.584, -3.751], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [1.372, 6.584, 0.538], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [2.616, 6.584, -1.607], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-4.77, 9.105, -2.526], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-4.48, 8.968, -0.796], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1))", pos = [-5.971, 9.73, -1.379], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 2))", pos = [-3.995, 6.507, -1.524], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3))", pos = [-3.131, 6.584, -1.606], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3))", pos = [-1.24, 6.584, -1.606], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3))", pos = [-5.93, 7.663, -1.524], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3))", pos = [-0.557, 6.584, -2.82], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3))", pos = [-0.557, 6.584, -0.392], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3))", pos = [0.836, 6.584, -2.811], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3))", pos = [0.836, 6.584, -0.401], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3))", pos = [1.535, 6.584, -1.606], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 4))", pos = [-5.222, 6.507, -1.525], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 4))", pos = [-6.028, 5.232, -1.518], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 4))", pos = [-5.249, 8.941, -1.563], state = 1)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [-7.558, 3.378, 0.726], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [-3.745, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [-10.493, 4.495, 0.675], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [-1.841, 1.309, 1.228], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [-1.841, 1.309, -3.06], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [0.636, 1.309, -3.061], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [1.881, 1.309, -0.917], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-9.534, 6.601, 0.871], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-8.043, 5.839, 1.454], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-8.333, 5.976, -0.276], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-10.664, 2.274, 0.762], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-9.295, 1.514, 1.596], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-9.342, 1.532, -0.159], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 6))", pos = [-1.292, 1.309, 0.298], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 6))", pos = [0.8, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 7))", pos = [-1.975, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 7))", pos = [-1.292, 1.309, -2.13], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 7))", pos = [0.101, 1.309, -2.121], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 8))", pos = [-8.785, 3.378, 0.725], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 8))", pos = [-9.493, 4.534, 0.726], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 4), (0, 6))", pos = [-9.591, 2.103, 0.732], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 4), (0, 6))", pos = [-8.812, 5.812, 0.687], state = 2)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [-6.93, 7.624, -1.575], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [-1.106, 6.584, -3.75], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [-1.106, 6.584, 0.538], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [1.371, 6.584, -3.751], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [1.372, 6.584, 0.538], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3))", pos = [2.616, 6.584, -1.607], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-5.779, 4.661, -2.409], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-5.732, 4.643, -0.654], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-7.101, 5.403, -1.488], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-4.77, 9.105, -2.526], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-4.48, 8.968, -0.796], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4))", pos = [-5.971, 9.73, -1.379], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7))", pos = [-0.557, 6.584, -2.82], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7))", pos = [-0.557, 6.584, -0.392], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7))", pos = [0.836, 6.584, -2.811], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7))", pos = [0.836, 6.584, -0.401], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7))", pos = [1.535, 6.584, -1.606], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 9))", pos = [-3.131, 6.584, -1.606], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 9))", pos = [-1.24, 6.584, -1.606], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 9))", pos = [-5.93, 7.663, -1.524], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 4), (0, 6))", pos = [-5.249, 8.941, -1.563], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 4), (0, 7))", pos = [-6.028, 5.232, -1.518], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 4), (0, 12))", pos = [-5.222, 6.507, -1.525], state = 2)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 6))", pos = [-1.841, 1.309, 1.228], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 6))", pos = [1.881, 1.309, -0.917], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7))", pos = [-3.745, 1.309, -0.916], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7))", pos = [-1.841, 1.309, -3.06], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7))", pos = [0.636, 1.309, -3.061], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 8))", pos = [-7.558, 3.378, 0.726], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 8))", pos = [-10.493, 4.495, 0.675], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-9.534, 6.601, 0.871], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-8.043, 5.839, 1.454], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-8.333, 5.976, -0.276], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-10.664, 2.274, 0.762], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-9.295, 1.514, 1.596], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-9.342, 1.532, -0.159], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 7), (0, 16))", pos = [-1.975, 1.309, -0.916], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 7), (0, 16))", pos = [0.101, 1.309, -2.121], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 7), (0, 17))", pos = [-1.292, 1.309, -2.13], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 4), (0, 6), (0, 20))", pos = [-9.591, 2.103, 0.732], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 4), (0, 6), (0, 20))", pos = [-8.812, 5.812, 0.687], state = 3)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7))", pos = [-1.106, 6.584, -3.75], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7))", pos = [-1.106, 6.584, 0.538], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7))", pos = [1.371, 6.584, -3.751], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7))", pos = [1.372, 6.584, 0.538], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7))", pos = [2.616, 6.584, -1.607], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 9))", pos = [-6.93, 7.624, -1.575], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-4.77, 9.105, -2.526], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-4.48, 8.968, -0.796], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6))", pos = [-5.971, 9.73, -1.379], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 7))", pos = [-5.779, 4.661, -2.409], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 7))", pos = [-5.732, 4.643, -0.654], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 7))", pos = [-7.101, 5.403, -1.488], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7), (0, 17))", pos = [0.836, 6.584, -2.811], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7), (0, 17))", pos = [0.836, 6.584, -0.401], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7), (0, 17))", pos = [1.535, 6.584, -1.606], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7), (0, 19))", pos = [-0.557, 6.584, -2.82], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7), (0, 19))", pos = [-0.557, 6.584, -0.392], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 4), (0, 6), (0, 21))", pos = [-5.249, 8.941, -1.563], state = 3)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 16))", pos = [-3.745, 1.309, -0.916], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 16))", pos = [0.636, 1.309, -3.061], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 17))", pos = [-1.841, 1.309, -3.06], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 20))", pos = [-9.534, 6.601, 0.871], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 20))", pos = [-8.043, 5.839, 1.454], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 20))", pos = [-8.333, 5.976, -0.276], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 20))", pos = [-10.664, 2.274, 0.762], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 20))", pos = [-9.295, 1.514, 1.596], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 20))", pos = [-9.342, 1.532, -0.159], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 3), (0, 7), (0, 17), (0, 39))", pos = [-1.292, 1.309, -2.13], state = 4)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 17))", pos = [1.371, 6.584, -3.751], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 17))", pos = [1.372, 6.584, 0.538], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 17))", pos = [2.616, 6.584, -1.607], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 19))", pos = [-1.106, 6.584, -3.75], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 19))", pos = [-1.106, 6.584, 0.538], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 21))", pos = [-4.77, 9.105, -2.526], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 21))", pos = [-4.48, 8.968, -0.796], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 4), (0, 6), (0, 21))", pos = [-5.971, 9.73, -1.379], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7), (0, 17), (0, 41))", pos = [1.535, 6.584, -1.606], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7), (0, 17), (0, 43))", pos = [0.836, 6.584, -2.811], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 3), (0, 7), (0, 17), (0, 43))", pos = [0.836, 6.584, -0.401], state = 4)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol1_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 17), (0, 39))", pos = [-1.841, 1.309, -3.06], state = 5)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 17), (0, 41))", pos = [2.616, 6.584, -1.607], state = 5)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 17), (0, 43))", pos = [1.371, 6.584, -3.751], state = 5)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol2_equiv", label="((0, 1), (0, 1), (0, 3), (0, 7), (0, 17), (0, 43))", pos = [1.372, 6.584, 0.538], state = 5)

cmd.set("label_size", 24, "mol2_equiv")

cmd.pseudoatom("mol1_equiv", label="[]", pos = [], state = 6)

cmd.set("label_size", 24, "mol1_equiv")

cmd.pseudoatom("mol2_equiv", label="[]", pos = [], state = 6)

cmd.set("label_size", 24, "mol2_equiv")


visualization_group = cmd.group("visualization_group")
cmd.group("visualization_group", "open")

cmd.group("visualization_group", "mol1_equiv", "add"),
cmd.group("visualization_group", "mol2_equiv", "add"),
cmd.group("visualization_group", "mol1_equiv", "add"),
cmd.group("visualization_group", "mol2_equiv", "add"),
cmd.group("visualization_group", "mol1_equiv", "add"),
cmd.group("visualization_group", "mol2_equiv", "add"),
cmd.group("visualization_group", "mol1_equiv", "add"),
cmd.group("visualization_group", "mol2_equiv", "add"),
cmd.group("visualization_group", "mol1_equiv", "add"),
cmd.group("visualization_group", "mol2_equiv", "add"),
cmd.group("visualization_group", "mol1_equiv", "add"),
cmd.group("visualization_group", "mol2_equiv", "add")

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

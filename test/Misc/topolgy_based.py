
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


mol1 = [
        
COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,0.16,1.309,0.392,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.975,1.309,-0.916,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.292,1.309,0.298,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.292,1.309,-2.13,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,0.101,1.309,-2.121,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,0.8,1.309,-0.916,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-3.745,1.309,-0.916,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.841,1.309,1.228,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.841,1.309,-3.06,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,0.636,1.309,-3.061,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.881,1.309,-0.917,0.3

            ]
cmd.load_cgo(mol1, "mol1", state=1)
cmd.set("cgo_transparency", 0, "mol1")
        
cmd.pseudoatom("mol1_labels", label="1.0", pos = [-8.785, 3.378, 0.725], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [0.16, 1.309, 0.392], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-9.493, 4.534, 0.726], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-7.558, 3.378, 0.726], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-9.591, 2.103, 0.732], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-8.812, 5.812, 0.687], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-1.975, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-1.292, 1.309, 0.298], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-1.292, 1.309, -2.13], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [0.101, 1.309, -2.121], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [0.8, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-3.745, 1.309, -0.916], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-10.493, 4.495, 0.675], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-9.534, 6.601, 0.871], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-8.043, 5.839, 1.454], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-8.333, 5.976, -0.276], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-10.664, 2.274, 0.762], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-9.295, 1.514, 1.596], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-9.342, 1.532, -0.159], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-1.841, 1.309, 1.228], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-1.841, 1.309, -3.06], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [0.636, 1.309, -3.061], state = 1)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [1.881, 1.309, -0.917], state = 1)

cmd.set("label_size", 24, "mol1_labels")


mol2 = [
        
COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-3.995,6.507,-1.524,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-3.131,6.584,-1.606,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.222,6.507,-1.525,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.24,6.584,-1.606,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.93,7.663,-1.524,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-0.557,6.584,-2.82,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-0.557,6.584,-0.392,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-6.93,7.624,-1.575,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.106,6.584,-3.75,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.106,6.584,0.538,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-6.028,5.232,-1.518,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.779,4.661,-2.409,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.732,4.643,-0.654,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-7.101,5.403,-1.488,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,0.836,6.584,-2.811,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,0.836,6.584,-0.401,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.371,6.584,-3.751,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.372,6.584,0.538,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.535,6.584,-1.606,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,2.616,6.584,-1.607,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.249,8.941,-1.563,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.77,9.105,-2.526,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.48,8.968,-0.796,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.971,9.73,-1.379,0.3

            ]
cmd.load_cgo(mol2, "mol2", state=1)
cmd.set("cgo_transparency", 0, "mol2")
        
cmd.pseudoatom("mol2_labels", label="1.0", pos = [-3.995, 6.507, -1.524], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-3.131, 6.584, -1.606], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.222, 6.507, -1.525], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-1.24, 6.584, -1.606], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.93, 7.663, -1.524], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-0.557, 6.584, -2.82], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-0.557, 6.584, -0.392], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-6.93, 7.624, -1.575], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-1.106, 6.584, -3.75], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-1.106, 6.584, 0.538], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-6.028, 5.232, -1.518], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.779, 4.661, -2.409], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.732, 4.643, -0.654], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-7.101, 5.403, -1.488], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [0.836, 6.584, -2.811], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [0.836, 6.584, -0.401], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [1.371, 6.584, -3.751], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [1.372, 6.584, 0.538], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [1.535, 6.584, -1.606], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [2.616, 6.584, -1.607], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.249, 8.941, -1.563], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-4.77, 9.105, -2.526], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-4.48, 8.968, -0.796], state = 1)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.971, 9.73, -1.379], state = 1)

cmd.set("label_size", 24, "mol2_labels")


mol1 = [
        
COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.7398692810457518,0.884967320261438,0.9333333333333333,1.0,SPHERE,0.16,1.309,0.392,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-1.975,1.309,-0.916,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-1.292,1.309,0.298,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-1.292,1.309,-2.13,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,0.101,1.309,-2.121,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,0.8,1.309,-0.916,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-3.745,1.309,-0.916,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.841,1.309,1.228,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.841,1.309,-3.06,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,0.636,1.309,-3.061,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.881,1.309,-0.917,0.3

            ]
cmd.load_cgo(mol1, "mol1", state=2)
cmd.set("cgo_transparency", 0, "mol1")
        
cmd.pseudoatom("mol1_labels", label="3.0", pos = [-8.785, 3.378, 0.725], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="2.0", pos = [0.16, 1.309, 0.392], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-9.493, 4.534, 0.726], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-7.558, 3.378, 0.726], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="4.0", pos = [-9.591, 2.103, 0.732], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="4.0", pos = [-8.812, 5.812, 0.687], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-1.975, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-1.292, 1.309, 0.298], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-1.292, 1.309, -2.13], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [0.101, 1.309, -2.121], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [0.8, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-3.745, 1.309, -0.916], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-10.493, 4.495, 0.675], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-9.534, 6.601, 0.871], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-8.043, 5.839, 1.454], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-8.333, 5.976, -0.276], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-10.664, 2.274, 0.762], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-9.295, 1.514, 1.596], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-9.342, 1.532, -0.159], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-1.841, 1.309, 1.228], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [-1.841, 1.309, -3.06], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [0.636, 1.309, -3.061], state = 2)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="1.0", pos = [1.881, 1.309, -0.917], state = 2)

cmd.set("label_size", 24, "mol1_labels")


mol2 = [
        
COLOR,0.7398692810457518,0.884967320261438,0.9333333333333333,1.0,SPHERE,-3.995,6.507,-1.524,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-3.131,6.584,-1.606,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-5.222,6.507,-1.525,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-1.24,6.584,-1.606,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-5.93,7.663,-1.524,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-0.557,6.584,-2.82,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-0.557,6.584,-0.392,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-6.93,7.624,-1.575,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.106,6.584,-3.75,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.106,6.584,0.538,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-6.028,5.232,-1.518,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.779,4.661,-2.409,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.732,4.643,-0.654,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-7.101,5.403,-1.488,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,0.836,6.584,-2.811,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,0.836,6.584,-0.401,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.371,6.584,-3.751,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.372,6.584,0.538,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,1.535,6.584,-1.606,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,2.616,6.584,-1.607,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-5.249,8.941,-1.563,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.77,9.105,-2.526,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.48,8.968,-0.796,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.971,9.73,-1.379,0.3

            ]
cmd.load_cgo(mol2, "mol2", state=2)
cmd.set("cgo_transparency", 0, "mol2")
        
cmd.pseudoatom("mol2_labels", label="2.0", pos = [-3.995, 6.507, -1.524], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [-3.131, 6.584, -1.606], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-5.222, 6.507, -1.525], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [-1.24, 6.584, -1.606], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [-5.93, 7.663, -1.524], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [-0.557, 6.584, -2.82], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [-0.557, 6.584, -0.392], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-6.93, 7.624, -1.575], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-1.106, 6.584, -3.75], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-1.106, 6.584, 0.538], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-6.028, 5.232, -1.518], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.779, 4.661, -2.409], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.732, 4.643, -0.654], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-7.101, 5.403, -1.488], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [0.836, 6.584, -2.811], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [0.836, 6.584, -0.401], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [1.371, 6.584, -3.751], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [1.372, 6.584, 0.538], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [1.535, 6.584, -1.606], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [2.616, 6.584, -1.607], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-5.249, 8.941, -1.563], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-4.77, 9.105, -2.526], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-4.48, 8.968, -0.796], state = 2)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="1.0", pos = [-5.971, 9.73, -1.379], state = 2)

cmd.set("label_size", 24, "mol2_labels")


mol1 = [
        
COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.996078431372549,0.8784313725490196,0.5647058823529412,1.0,SPHERE,0.16,1.309,0.392,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.996078431372549,0.8784313725490196,0.5647058823529412,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.996078431372549,0.8784313725490196,0.5647058823529412,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.9568627450980393,0.42745098039215684,0.2627450980392157,1.0,SPHERE,-1.975,1.309,-0.916,0.3,COLOR,0.996078431372549,0.8784313725490196,0.5647058823529412,1.0,SPHERE,-1.292,1.309,0.298,0.3,COLOR,0.9568627450980393,0.42745098039215684,0.2627450980392157,1.0,SPHERE,-1.292,1.309,-2.13,0.3,COLOR,0.9568627450980393,0.42745098039215684,0.2627450980392157,1.0,SPHERE,0.101,1.309,-2.121,0.3,COLOR,0.996078431372549,0.8784313725490196,0.5647058823529412,1.0,SPHERE,0.8,1.309,-0.916,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-3.745,1.309,-0.916,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.45490196078431383,0.6784313725490198,0.8196078431372549,1.0,SPHERE,-9.534,6.601,0.871,0.3,COLOR,0.45490196078431383,0.6784313725490198,0.8196078431372549,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.45490196078431383,0.6784313725490198,0.8196078431372549,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.45490196078431383,0.6784313725490198,0.8196078431372549,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.45490196078431383,0.6784313725490198,0.8196078431372549,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.45490196078431383,0.6784313725490198,0.8196078431372549,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.841,1.309,1.228,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.841,1.309,-3.06,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,0.636,1.309,-3.061,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.881,1.309,-0.917,0.3

            ]
cmd.load_cgo(mol1, "mol1", state=3)
cmd.set("cgo_transparency", 0, "mol1")
        
cmd.pseudoatom("mol1_labels", label="8.0", pos = [-8.785, 3.378, 0.725], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [0.16, 1.309, 0.392], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="8.0", pos = [-9.493, 4.534, 0.726], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-7.558, 3.378, 0.726], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-9.591, 2.103, 0.732], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-8.812, 5.812, 0.687], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="7.0", pos = [-1.975, 1.309, -0.916], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-1.292, 1.309, 0.298], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="7.0", pos = [-1.292, 1.309, -2.13], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="7.0", pos = [0.101, 1.309, -2.121], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [0.8, 1.309, -0.916], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-3.745, 1.309, -0.916], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-10.493, 4.495, 0.675], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="4.0", pos = [-9.534, 6.601, 0.871], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="4.0", pos = [-8.043, 5.839, 1.454], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="4.0", pos = [-8.333, 5.976, -0.276], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="4.0", pos = [-10.664, 2.274, 0.762], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="4.0", pos = [-9.295, 1.514, 1.596], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="4.0", pos = [-9.342, 1.532, -0.159], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-1.841, 1.309, 1.228], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [-1.841, 1.309, -3.06], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [0.636, 1.309, -3.061], state = 3)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="3.0", pos = [1.881, 1.309, -0.917], state = 3)

cmd.set("label_size", 24, "mol1_labels")


mol2 = [
        
COLOR,0.9308727412533642,0.9732410611303345,0.8761245674740483,1.0,SPHERE,-3.995,6.507,-1.524,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-3.131,6.584,-1.606,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-5.222,6.507,-1.525,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-1.24,6.584,-1.606,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-5.93,7.663,-1.524,0.3,COLOR,0.9308727412533642,0.9732410611303345,0.8761245674740483,1.0,SPHERE,-0.557,6.584,-2.82,0.3,COLOR,0.9308727412533642,0.9732410611303345,0.8761245674740483,1.0,SPHERE,-0.557,6.584,-0.392,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-6.93,7.624,-1.575,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.106,6.584,-3.75,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.106,6.584,0.538,0.3,COLOR,0.9308727412533642,0.9732410611303345,0.8761245674740483,1.0,SPHERE,-6.028,5.232,-1.518,0.3,COLOR,0.28865820838139183,0.48035371011149564,0.7170319108035371,1.0,SPHERE,-5.779,4.661,-2.409,0.3,COLOR,0.28865820838139183,0.48035371011149564,0.7170319108035371,1.0,SPHERE,-5.732,4.643,-0.654,0.3,COLOR,0.28865820838139183,0.48035371011149564,0.7170319108035371,1.0,SPHERE,-7.101,5.403,-1.488,0.3,COLOR,0.9308727412533642,0.9732410611303345,0.8761245674740483,1.0,SPHERE,0.836,6.584,-2.811,0.3,COLOR,0.9308727412533642,0.9732410611303345,0.8761245674740483,1.0,SPHERE,0.836,6.584,-0.401,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.371,6.584,-3.751,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.372,6.584,0.538,0.3,COLOR,0.9308727412533642,0.9732410611303345,0.8761245674740483,1.0,SPHERE,1.535,6.584,-1.606,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,2.616,6.584,-1.607,0.3,COLOR,0.7398692810457518,0.884967320261438,0.9333333333333333,1.0,SPHERE,-5.249,8.941,-1.563,0.3,COLOR,0.28865820838139183,0.48035371011149564,0.7170319108035371,1.0,SPHERE,-4.77,9.105,-2.526,0.3,COLOR,0.28865820838139183,0.48035371011149564,0.7170319108035371,1.0,SPHERE,-4.48,8.968,-0.796,0.3,COLOR,0.28865820838139183,0.48035371011149564,0.7170319108035371,1.0,SPHERE,-5.971,9.73,-1.379,0.3

            ]
cmd.load_cgo(mol2, "mol2", state=3)
cmd.set("cgo_transparency", 0, "mol2")
        
cmd.pseudoatom("mol2_labels", label="7.0", pos = [-3.995, 6.507, -1.524], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="9.0", pos = [-3.131, 6.584, -1.606], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="12.0", pos = [-5.222, 6.507, -1.525], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="9.0", pos = [-1.24, 6.584, -1.606], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="9.0", pos = [-5.93, 7.663, -1.524], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [-0.557, 6.584, -2.82], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [-0.557, 6.584, -0.392], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [-6.93, 7.624, -1.575], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [-1.106, 6.584, -3.75], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [-1.106, 6.584, 0.538], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [-6.028, 5.232, -1.518], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-5.779, 4.661, -2.409], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-5.732, 4.643, -0.654], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-7.101, 5.403, -1.488], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [0.836, 6.584, -2.811], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [0.836, 6.584, -0.401], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [1.371, 6.584, -3.751], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [1.372, 6.584, 0.538], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [1.535, 6.584, -1.606], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="3.0", pos = [2.616, 6.584, -1.607], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="6.0", pos = [-5.249, 8.941, -1.563], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-4.77, 9.105, -2.526], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-4.48, 8.968, -0.796], state = 3)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="4.0", pos = [-5.971, 9.73, -1.379], state = 3)

cmd.set("label_size", 24, "mol2_labels")


mol1 = [
        
COLOR,0.9610149942329873,0.457439446366782,0.2765859284890427,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.9118031526336026,0.9658592848904267,0.9111880046136099,1.0,SPHERE,0.16,1.309,0.392,0.3,COLOR,0.9610149942329873,0.457439446366782,0.2765859284890427,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.34648212226066905,0.5492502883506345,0.7527104959630911,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.9873125720876587,0.647366397539408,0.36424452133794694,1.0,SPHERE,-1.975,1.309,-0.916,0.3,COLOR,0.9873125720876587,0.647366397539408,0.36424452133794694,1.0,SPHERE,-1.292,1.309,0.298,0.3,COLOR,0.9610149942329873,0.457439446366782,0.2765859284890427,1.0,SPHERE,-1.292,1.309,-2.13,0.3,COLOR,0.9873125720876587,0.647366397539408,0.36424452133794694,1.0,SPHERE,0.101,1.309,-2.121,0.3,COLOR,0.9873125720876587,0.647366397539408,0.36424452133794694,1.0,SPHERE,0.8,1.309,-0.916,0.3,COLOR,0.247520184544406,0.38615916955017304,0.6701268742791234,1.0,SPHERE,-3.745,1.309,-0.916,0.3,COLOR,0.34648212226066905,0.5492502883506345,0.7527104959630911,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-1.841,1.309,1.228,0.3,COLOR,0.247520184544406,0.38615916955017304,0.6701268742791234,1.0,SPHERE,-1.841,1.309,-3.06,0.3,COLOR,0.247520184544406,0.38615916955017304,0.6701268742791234,1.0,SPHERE,0.636,1.309,-3.061,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,1.881,1.309,-0.917,0.3

            ]
cmd.load_cgo(mol1, "mol1", state=4)
cmd.set("cgo_transparency", 0, "mol1")
        
cmd.pseudoatom("mol1_labels", label="17.0", pos = [-8.785, 3.378, 0.725], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="12.0", pos = [0.16, 1.309, 0.392], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="17.0", pos = [-9.493, 4.534, 0.726], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="8.0", pos = [-7.558, 3.378, 0.726], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="20.0", pos = [-9.591, 2.103, 0.732], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="20.0", pos = [-8.812, 5.812, 0.687], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="16.0", pos = [-1.975, 1.309, -0.916], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="16.0", pos = [-1.292, 1.309, 0.298], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="17.0", pos = [-1.292, 1.309, -2.13], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="16.0", pos = [0.101, 1.309, -2.121], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="16.0", pos = [0.8, 1.309, -0.916], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="7.0", pos = [-3.745, 1.309, -0.916], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="8.0", pos = [-10.493, 4.495, 0.675], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-9.534, 6.601, 0.871], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-8.043, 5.839, 1.454], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-8.333, 5.976, -0.276], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-10.664, 2.274, 0.762], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-9.295, 1.514, 1.596], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-9.342, 1.532, -0.159], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [-1.841, 1.309, 1.228], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="7.0", pos = [-1.841, 1.309, -3.06], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="7.0", pos = [0.636, 1.309, -3.061], state = 4)

cmd.set("label_size", 24, "mol1_labels")

cmd.pseudoatom("mol1_labels", label="6.0", pos = [1.881, 1.309, -0.917], state = 4)

cmd.set("label_size", 24, "mol1_labels")


mol2 = [
        
COLOR,0.9970011534025375,0.9070357554786621,0.6080738177623991,1.0,SPHERE,-3.995,6.507,-1.524,0.3,COLOR,0.9033448673587082,0.314878892733564,0.2110726643598616,1.0,SPHERE,-3.131,6.584,-1.606,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-5.222,6.507,-1.525,0.3,COLOR,0.993925413302576,0.7707804690503652,0.463514033064206,1.0,SPHERE,-1.24,6.584,-1.606,0.3,COLOR,0.9970011534025375,0.9070357554786621,0.6080738177623991,1.0,SPHERE,-5.93,7.663,-1.524,0.3,COLOR,0.9999231064975009,0.9976163014225298,0.7454056132256824,1.0,SPHERE,-0.557,6.584,-2.82,0.3,COLOR,0.9999231064975009,0.9976163014225298,0.7454056132256824,1.0,SPHERE,-0.557,6.584,-0.392,0.3,COLOR,0.29588619761630147,0.488965782391388,0.7214917339484814,1.0,SPHERE,-6.93,7.624,-1.575,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,-1.106,6.584,-3.75,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,-1.106,6.584,0.538,0.3,COLOR,0.9923875432525952,0.6938869665513264,0.39123414071510954,1.0,SPHERE,-6.028,5.232,-1.518,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,-5.779,4.661,-2.409,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,-5.732,4.643,-0.654,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,-7.101,5.403,-1.488,0.3,COLOR,0.9070357554786621,0.9640138408304498,0.9199538638985004,1.0,SPHERE,0.836,6.584,-2.811,0.3,COLOR,0.9070357554786621,0.9640138408304498,0.9199538638985004,1.0,SPHERE,0.836,6.584,-0.401,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,1.371,6.584,-3.751,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,1.372,6.584,0.538,0.3,COLOR,0.9070357554786621,0.9640138408304498,0.9199538638985004,1.0,SPHERE,1.535,6.584,-1.606,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,2.616,6.584,-1.607,0.3,COLOR,0.9970011534025375,0.9070357554786621,0.6080738177623991,1.0,SPHERE,-5.249,8.941,-1.563,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.77,9.105,-2.526,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.48,8.968,-0.796,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-5.971,9.73,-1.379,0.3

            ]
cmd.load_cgo(mol2, "mol2", state=4)
cmd.set("cgo_transparency", 0, "mol2")
        
cmd.pseudoatom("mol2_labels", label="21.0", pos = [-3.995, 6.507, -1.524], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="28.0", pos = [-3.131, 6.584, -1.606], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="32.0", pos = [-5.222, 6.507, -1.525], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="23.0", pos = [-1.24, 6.584, -1.606], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="21.0", pos = [-5.93, 7.663, -1.524], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="19.0", pos = [-0.557, 6.584, -2.82], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="19.0", pos = [-0.557, 6.584, -0.392], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="9.0", pos = [-6.93, 7.624, -1.575], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [-1.106, 6.584, -3.75], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [-1.106, 6.584, 0.538], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="24.0", pos = [-6.028, 5.232, -1.518], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [-5.779, 4.661, -2.409], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [-5.732, 4.643, -0.654], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [-7.101, 5.403, -1.488], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="17.0", pos = [0.836, 6.584, -2.811], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="17.0", pos = [0.836, 6.584, -0.401], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [1.371, 6.584, -3.751], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [1.372, 6.584, 0.538], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="17.0", pos = [1.535, 6.584, -1.606], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="7.0", pos = [2.616, 6.584, -1.607], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="21.0", pos = [-5.249, 8.941, -1.563], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="6.0", pos = [-4.77, 9.105, -2.526], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="6.0", pos = [-4.48, 8.968, -0.796], state = 4)

cmd.set("label_size", 24, "mol2_labels")

cmd.pseudoatom("mol2_labels", label="6.0", pos = [-5.971, 9.73, -1.379], state = 4)

cmd.set("label_size", 24, "mol2_labels")


visualization_group = cmd.group("visualization_group")
cmd.group("visualization_group", "open")

cmd.group("visualization_group", "mol1", "add"),
cmd.group("visualization_group", "mol1_labels", "add"),
cmd.group("visualization_group", "mol2", "add"),
cmd.group("visualization_group", "mol2_labels", "add"),
cmd.group("visualization_group", "mol1", "add"),
cmd.group("visualization_group", "mol1_labels", "add"),
cmd.group("visualization_group", "mol2", "add"),
cmd.group("visualization_group", "mol2_labels", "add"),
cmd.group("visualization_group", "mol1", "add"),
cmd.group("visualization_group", "mol1_labels", "add"),
cmd.group("visualization_group", "mol2", "add"),
cmd.group("visualization_group", "mol2_labels", "add"),
cmd.group("visualization_group", "mol1", "add"),
cmd.group("visualization_group", "mol1_labels", "add"),
cmd.group("visualization_group", "mol2", "add"),
cmd.group("visualization_group", "mol2_labels", "add")

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()


from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


history_1 = [
        
COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=1)
cmd.set("cgo_transparency", 0, "history_1")
        

history_1 = [
        
COLOR,0.4633602460592081,0.6851980007689351,0.8232987312572088,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.41153402537485584,0.6267589388696656,0.7928489042675894,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.3537101114955786,0.5578623606305267,0.7571703191080355,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=2)
cmd.set("cgo_transparency", 0, "history_1")
        

history_1 = [
        
COLOR,0.9953094963475586,0.8399846212995001,0.5285659361783929,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.8295271049596311,0.9289504036908882,0.9587081891580161,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.9070357554786621,0.9640138408304498,0.9199538638985004,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.31034217608612075,0.5061899269511727,0.7304113802383699,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.2659746251441753,0.4442906574394464,0.6987312572087659,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.2659746251441753,0.4442906574394464,0.6987312572087659,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.21368704344482892,0.2795847750865052,0.617685505574779,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.2659746251441753,0.4442906574394464,0.6987312572087659,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.2659746251441753,0.4442906574394464,0.6987312572087659,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.2659746251441753,0.4442906574394464,0.6987312572087659,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.22599000384467513,0.31833910034602075,0.6367550941945406,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=3)
cmd.set("cgo_transparency", 0, "history_1")
        

history_1 = [
        
COLOR,0.8587466359092657,0.22106881968473663,0.16801230296039987,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.9982314494425221,0.9451749327181853,0.6658977316416763,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.9982314494425221,0.9451749327181853,0.6658977316416763,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.3898500576701269,0.6009227220299886,0.7794694348327567,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.38262206843521723,0.5923106497500962,0.7750096116878125,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.38262206843521723,0.5923106497500962,0.7750096116878125,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.2536716647443291,0.4055363321799308,0.6796616685890043,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.22291426374471357,0.3086505190311419,0.6319876970396002,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.22291426374471357,0.3086505190311419,0.6319876970396002,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.2413687043444829,0.3667820069204153,0.6605920799692426,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.25982314494425224,0.42491349480968865,0.6891964628988851,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.25982314494425224,0.42491349480968865,0.6891964628988851,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.22291426374471357,0.3086505190311419,0.6319876970396002,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.22291426374471357,0.3086505190311419,0.6319876970396002,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.25982314494425224,0.42491349480968865,0.6891964628988851,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.22291426374471357,0.3086505190311419,0.6319876970396002,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.21061130334486736,0.2698961937716263,0.6129181084198385,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.1952326028450596,0.22145328719723184,0.5890811226451366,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.1952326028450596,0.22145328719723184,0.5890811226451366,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.1952326028450596,0.22145328719723184,0.5890811226451366,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=4)
cmd.set("cgo_transparency", 0, "history_1")
        

history_1 = [
        
COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.9943867743175702,0.7938485198000771,0.4851980007689352,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.9965397923875433,0.8927335640138409,0.5863898500576701,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.4887351018838909,0.7054978854286814,0.8343713956170704,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.4887351018838909,0.7054978854286814,0.8343713956170704,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.4887351018838909,0.7054978854286814,0.8343713956170704,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.32479815455594,0.5234140715109573,0.7393310265282584,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.28143021914648214,0.4717416378316033,0.7125720876585929,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.28143021914648214,0.4717416378316033,0.7125720876585929,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.2742022299115725,0.4631295655517109,0.7081122645136487,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.26289888504421377,0.4346020761245675,0.6939638600538255,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.26289888504421377,0.4346020761245675,0.6939638600538255,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.20753556324490582,0.2602076124567474,0.6081507112648982,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.20753556324490582,0.2602076124567474,0.6081507112648982,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.20753556324490582,0.2602076124567474,0.6081507112648982,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.23829296424452134,0.35709342560553636,0.6558246828143023,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.23829296424452134,0.35709342560553636,0.6558246828143023,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.25059592464436753,0.395847750865052,0.6748942714340639,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.23829296424452134,0.35709342560553636,0.6558246828143023,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.21061130334486736,0.2698961937716263,0.6129181084198385,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=5)
cmd.set("cgo_transparency", 0, "history_1")
        

history_1 = [
        
COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.7316416762783547,0.08119953863898521,0.15071126489811612,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.9934640522875817,0.7477124183006536,0.4418300653594771,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.9968473663975395,0.9022683583237218,0.6008458285274896,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.5733179546328336,0.7731641676278356,0.871280276816609,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.5648596693579393,0.7663975394079201,0.8675893886966551,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.5648596693579393,0.7663975394079201,0.8675893886966551,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.4187620146097656,0.635371011149558,0.7973087274125337,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.3898500576701269,0.6009227220299886,0.7794694348327567,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.3898500576701269,0.6009227220299886,0.7794694348327567,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.34648212226066905,0.5492502883506345,0.7527104959630911,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.2742022299115725,0.4631295655517109,0.7081122645136487,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.2742022299115725,0.4631295655517109,0.7081122645136487,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.23521722414455978,0.34740484429065743,0.6510572856593618,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.23521722414455978,0.34740484429065743,0.6510572856593618,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.23521722414455978,0.34740484429065743,0.6510572856593618,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.247520184544406,0.38615916955017304,0.6701268742791234,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.247520184544406,0.38615916955017304,0.6701268742791234,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.25059592464436753,0.395847750865052,0.6748942714340639,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.23829296424452134,0.35709342560553636,0.6558246828143023,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.21368704344482892,0.2795847750865052,0.617685505574779,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=6)
cmd.set("cgo_transparency", 0, "history_1")
        

history_1 = [
        
COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.7623990772779701,0.11072664359861592,0.15132641291810842,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.9925413302575933,0.7015763168012303,0.3984621299500192,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.9966935793925413,0.8975009611687812,0.5936178392925797,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.6494425221068819,0.8340638216070742,0.9044982698961938,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.6325259515570935,0.8205305651672434,0.8971164936562861,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.6325259515570935,0.8205305651672434,0.8971164936562861,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.5141099577085737,0.7257977700884276,0.845444059976932,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.4971933871587852,0.7122645136485968,0.8380622837370243,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.4971933871587852,0.7122645136485968,0.8380622837370243,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.4187620146097656,0.635371011149558,0.7973087274125337,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.29588619761630147,0.488965782391388,0.7214917339484814,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.29588619761630147,0.488965782391388,0.7214917339484814,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.2659746251441753,0.4442906574394464,0.6987312572087659,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.2659746251441753,0.4442906574394464,0.6987312572087659,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.2659746251441753,0.4442906574394464,0.6987312572087659,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.2536716647443291,0.4055363321799308,0.6796616685890043,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.2536716647443291,0.4055363321799308,0.6796616685890043,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.25059592464436753,0.395847750865052,0.6748942714340639,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.23829296424452134,0.35709342560553636,0.6558246828143023,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.21983852364475204,0.298961937716263,0.6272202998846598,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=7)
cmd.set("cgo_transparency", 0, "history_1")
        

history_1 = [
        
COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.7700884275278739,0.11810841983852365,0.15148019992310652,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.9900807381776241,0.6673587081891583,0.3734717416378317,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.9963860053825452,0.8879661668589004,0.5791618608227604,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.7154171472510573,0.8729719338715878,0.9264129181084199,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.690965013456363,0.8609765474817378,0.9194925028835064,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.690965013456363,0.8609765474817378,0.9194925028835064,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.6071510957324107,0.8002306805074972,0.8860438292964244,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.5902345251826222,0.7866974240676664,0.8786620530565168,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.5902345251826222,0.7866974240676664,0.8786620530565168,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.4887351018838909,0.7054978854286814,0.8343713956170704,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.31757016532103044,0.5148019992310651,0.7348712033833141,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.31757016532103044,0.5148019992310651,0.7348712033833141,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.3320261437908497,0.5320261437908498,0.7437908496732026,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.3320261437908497,0.5320261437908498,0.7437908496732026,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.3320261437908497,0.5320261437908498,0.7437908496732026,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.26289888504421377,0.4346020761245675,0.6939638600538255,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.26289888504421377,0.4346020761245675,0.6939638600538255,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.25059592464436753,0.395847750865052,0.6748942714340639,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.23521722414455978,0.34740484429065743,0.6510572856593618,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.22291426374471357,0.3086505190311419,0.6319876970396002,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=8)
cmd.set("cgo_transparency", 0, "history_1")
        

history_1 = [
        
COLOR,0.6470588235294118,0.0,0.14901960784313725,1.0,SPHERE,-7.558,3.378,0.726,0.3,COLOR,0.7777777777777778,0.12549019607843137,0.1516339869281046,1.0,SPHERE,-6.694281,3.455027,0.64431,0.3,COLOR,0.9845444059976932,0.6273740868896578,0.3550173010380623,1.0,SPHERE,-8.785,3.378,0.725,0.3,COLOR,0.996078431372549,0.8784313725490196,0.5647058823529412,1.0,SPHERE,-4.802274,3.455027,0.64431,0.3,COLOR,0.7724721261053442,0.900961168781238,0.942560553633218,1.0,SPHERE,-9.493,4.534,0.726,0.3,COLOR,0.7480199923106499,0.888965782391388,0.9356401384083045,1.0,SPHERE,-4.119275,3.455027,-0.569685,0.3,COLOR,0.7480199923106499,0.888965782391388,0.9356401384083045,1.0,SPHERE,-4.119275,3.455027,1.858314,0.3,COLOR,0.690965013456363,0.8609765474817378,0.9194925028835064,1.0,SPHERE,-10.493,4.495,0.675,0.3,COLOR,0.6663590926566706,0.8475970780469051,0.9118800461361015,1.0,SPHERE,-4.668274,3.455027,-1.499688,0.3,COLOR,0.6663590926566706,0.8475970780469051,0.9118800461361015,1.0,SPHERE,-4.668274,3.455027,2.78831,0.3,COLOR,0.5733179546328336,0.7731641676278356,0.871280276816609,1.0,SPHERE,-9.591,2.103,0.732,0.3,COLOR,0.33925413302575935,0.5406382160707421,0.748250672818147,1.0,SPHERE,-2.726279,3.455027,-0.560685,0.3,COLOR,0.33925413302575935,0.5406382160707421,0.748250672818147,1.0,SPHERE,-2.726279,3.455027,1.849314,0.3,COLOR,0.40430603613994626,0.6181468665897734,0.7883890811226452,1.0,SPHERE,-9.342,1.532,-0.159,0.3,COLOR,0.40430603613994626,0.6181468665897734,0.7883890811226452,1.0,SPHERE,-9.295,1.514,1.596,0.3,COLOR,0.40430603613994626,0.6181468665897734,0.7883890811226452,1.0,SPHERE,-10.664,2.274,0.762,0.3,COLOR,0.2690503652441369,0.45397923875432533,0.7034986543637063,1.0,SPHERE,-2.191279,3.455027,-1.500688,0.3,COLOR,0.2690503652441369,0.45397923875432533,0.7034986543637063,1.0,SPHERE,-2.190279,3.455027,2.78831,0.3,COLOR,0.2536716647443291,0.4055363321799308,0.6796616685890043,1.0,SPHERE,-2.027281,3.455027,0.64431,0.3,COLOR,0.23214148404459825,0.3377162629757786,0.6462898885044215,1.0,SPHERE,-0.946282,3.455027,0.64331,0.3,COLOR,0.22599000384467513,0.31833910034602075,0.6367550941945406,1.0,SPHERE,-8.812,5.812,0.687,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.333,5.976,-0.276,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-8.043,5.839,1.454,0.3,COLOR,0.19215686274509805,0.21176470588235294,0.5843137254901961,1.0,SPHERE,-9.534,6.601,0.871,0.3

            ]
cmd.load_cgo(history_1, "history_1", state=9)
cmd.set("cgo_transparency", 0, "history_1")
        

Group_0 = cmd.group("Group_0")
cmd.group("Group_0", "open")

cmd.group("Group_0", "history_1", "add"),
cmd.group("Group_0", "history_1", "add"),
cmd.group("Group_0", "history_1", "add"),
cmd.group("Group_0", "history_1", "add"),
cmd.group("Group_0", "history_1", "add"),
cmd.group("Group_0", "history_1", "add"),
cmd.group("Group_0", "history_1", "add"),
cmd.group("Group_0", "history_1", "add"),
cmd.group("Group_0", "history_1", "add")

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()
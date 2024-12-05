
from pymol.cgo import *
from pymol import cmd
import numpy as np
from chempy.brick import Brick
from collections import defaultdict
positions_viewport_callbacks = defaultdict(lambda: defaultdict(lambda: ViewportCallback([],0,0)))


cmd.set_color("dummy_0", (0.20445982314494426, 0.2505190311418685, 0.6033833141099577))
cmd.color("dummy_0", "test2 and state 1", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1")

cmd.set_color("dummy_1", (0.39707804690503656, 0.6095347943098809, 0.783929257977701))
cmd.color("dummy_1", "test2 and state 2", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 2")

cmd.set_color("dummy_2", (0.33925413302575935, 0.5406382160707421, 0.748250672818147))
cmd.color("dummy_2", "test2 and state 3", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 3")

cmd.set_color("dummy_3", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_3", "test2 and state 4", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 4")

cmd.set_color("dummy_4", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_4", "test2 and state 5", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 5")

cmd.set_color("dummy_5", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_5", "test2 and state 6", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 6")

cmd.set_color("dummy_6", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_6", "test2 and state 7", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 7")

cmd.set_color("dummy_7", (0.9033448673587082, 0.314878892733564, 0.2110726643598616))
cmd.color("dummy_7", "test2 and state 8", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 8")

cmd.set_color("dummy_8", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_8", "test2 and state 9", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 9")

cmd.set_color("dummy_9", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_9", "test2 and state 10", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 10")

cmd.set_color("dummy_10", (0.6828143021914649, 0.8569780853517878, 0.9171856978085352))
cmd.color("dummy_10", "test2 and state 11", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 11")

cmd.set_color("dummy_11", (0.9568627450980393, 0.42745098039215684, 0.2627450980392157))
cmd.color("dummy_11", "test2 and state 12", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 12")

cmd.set_color("dummy_12", (0.2659746251441753, 0.4442906574394464, 0.6987312572087659))
cmd.color("dummy_12", "test2 and state 13", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 13")

cmd.set_color("dummy_13", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_13", "test2 and state 14", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 14")

cmd.set_color("dummy_14", (0.5056516724336794, 0.7190311418685121, 0.8417531718569781))
cmd.color("dummy_14", "test2 and state 15", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 15")

cmd.set_color("dummy_15", (0.9256439830834294, 0.3617839292579777, 0.23260284505959247))
cmd.color("dummy_15", "test2 and state 16", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 16")

cmd.set_color("dummy_16", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_16", "test2 and state 17", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 17")

cmd.set_color("dummy_17", (0.22291426374471357, 0.3086505190311419, 0.6319876970396002))
cmd.color("dummy_17", "test2 and state 18", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 18")

cmd.set_color("dummy_18", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_18", "test2 and state 19", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 19")

cmd.set_color("dummy_19", (0.368166089965398, 0.5750865051903116, 0.766089965397924))
cmd.color("dummy_19", "test2 and state 20", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 20")

cmd.set_color("dummy_20", (0.9933102652825836, 0.7400230680507497, 0.43460207612456747))
cmd.color("dummy_20", "test2 and state 21", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 21")

cmd.set_color("dummy_21", (0.9953094963475586, 0.8399846212995001, 0.5285659361783929))
cmd.color("dummy_21", "test2 and state 22", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 22")

cmd.set_color("dummy_22", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_22", "test2 and state 23", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 23")

cmd.set_color("dummy_23", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_23", "test2 and state 24", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 24")

cmd.set_color("dummy_24", (0.38262206843521723, 0.5923106497500962, 0.7750096116878125))
cmd.color("dummy_24", "test2 and state 25", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 25")

cmd.set_color("dummy_25", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_25", "test2 and state 26", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 26")

cmd.set_color("dummy_26", (0.3031141868512111, 0.49757785467128035, 0.7259515570934256))
cmd.color("dummy_26", "test2 and state 27", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 27")

cmd.set_color("dummy_27", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_27", "test2 and state 28", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 28")

cmd.set_color("dummy_28", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_28", "test2 and state 29", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 29")

cmd.set_color("dummy_29", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_29", "test2 and state 30", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 30")

cmd.set_color("dummy_30", (0.9968473663975395, 0.9022683583237218, 0.6008458285274896))
cmd.color("dummy_30", "test2 and state 31", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 31")

cmd.set_color("dummy_31", (0.5479430988081508, 0.7528642829680893, 0.8602076124567475))
cmd.color("dummy_31", "test2 and state 32", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 32")

cmd.set_color("dummy_32", (0.9976163014225298, 0.9990772779700116, 0.7534025374855824))
cmd.color("dummy_32", "test2 and state 33", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 33")

cmd.set_color("dummy_33", (0.23829296424452134, 0.35709342560553636, 0.6558246828143023))
cmd.color("dummy_33", "test2 and state 34", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 34")

cmd.set_color("dummy_34", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_34", "test2 and state 35", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 35")

cmd.set_color("dummy_35", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_35", "test2 and state 36", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 36")

cmd.set_color("dummy_36", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_36", "test2 and state 37", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 37")

cmd.set_color("dummy_37", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_37", "test2 and state 38", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 38")

cmd.set_color("dummy_38", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_38", "test2 and state 39", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 39")

cmd.set_color("dummy_39", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_39", "test2 and state 40", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 40")

cmd.set_color("dummy_40", (0.9928489042675894, 0.7169550173010381, 0.4129181084198385))
cmd.color("dummy_40", "test2 and state 41", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 41")

cmd.set_color("dummy_41", (0.9873125720876587, 0.647366397539408, 0.36424452133794694))
cmd.color("dummy_41", "test2 and state 42", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 42")

cmd.set_color("dummy_42", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_42", "test2 and state 43", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 43")

cmd.set_color("dummy_43", (0.9873125720876587, 0.647366397539408, 0.36424452133794694))
cmd.color("dummy_43", "test2 and state 44", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 44")

cmd.set_color("dummy_44", (0.9968473663975395, 0.9022683583237218, 0.6008458285274896))
cmd.color("dummy_44", "test2 and state 45", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 45")

cmd.set_color("dummy_45", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_45", "test2 and state 46", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 46")

cmd.set_color("dummy_46", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_46", "test2 and state 47", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 47")

cmd.set_color("dummy_47", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_47", "test2 and state 48", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 48")

cmd.set_color("dummy_48", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_48", "test2 and state 49", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 49")

cmd.set_color("dummy_49", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_49", "test2 and state 50", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 50")

cmd.set_color("dummy_50", (0.9256439830834294, 0.3617839292579777, 0.23260284505959247))
cmd.color("dummy_50", "test2 and state 51", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 51")

cmd.set_color("dummy_51", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_51", "test2 and state 52", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 52")

cmd.set_color("dummy_52", (0.9985390234525182, 0.9547097270280661, 0.6803537101114956))
cmd.color("dummy_52", "test2 and state 53", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 53")

cmd.set_color("dummy_53", (0.9982314494425221, 0.9451749327181853, 0.6658977316416763))
cmd.color("dummy_53", "test2 and state 54", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 54")

cmd.set_color("dummy_54", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_54", "test2 and state 55", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 55")

cmd.set_color("dummy_55", (0.9951557093425606, 0.8322952710495963, 0.5213379469434832))
cmd.color("dummy_55", "test2 and state 56", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 56")

cmd.set_color("dummy_56", (0.6828143021914649, 0.8569780853517878, 0.9171856978085352))
cmd.color("dummy_56", "test2 and state 57", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 57")

cmd.set_color("dummy_57", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_57", "test2 and state 58", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 58")

cmd.set_color("dummy_58", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_58", "test2 and state 59", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 59")

cmd.set_color("dummy_59", (0.9997693194925029, 0.9928489042675894, 0.7381776239907728))
cmd.color("dummy_59", "test2 and state 60", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 60")

cmd.set_color("dummy_60", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_60", "test2 and state 61", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 61")

cmd.set_color("dummy_61", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_61", "test2 and state 62", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 62")

cmd.set_color("dummy_62", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_62", "test2 and state 63", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 63")

cmd.set_color("dummy_63", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_63", "test2 and state 64", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 64")

cmd.set_color("dummy_64", (0.9976163014225298, 0.9990772779700116, 0.7534025374855824))
cmd.color("dummy_64", "test2 and state 65", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 65")

cmd.set_color("dummy_65", (0.25674740484429065, 0.4152249134948097, 0.6844290657439447))
cmd.color("dummy_65", "test2 and state 66", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 66")

cmd.set_color("dummy_66", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_66", "test2 and state 67", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 67")

cmd.set_color("dummy_67", (0.4887351018838909, 0.7054978854286814, 0.8343713956170704))
cmd.color("dummy_67", "test2 and state 68", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 68")

cmd.set_color("dummy_68", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_68", "test2 and state 69", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 69")

cmd.set_color("dummy_69", (0.5986928104575164, 0.7934640522875818, 0.8823529411764706))
cmd.color("dummy_69", "test2 and state 70", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 70")

cmd.set_color("dummy_70", (0.37539407920030765, 0.5836985774702039, 0.7705497885428682))
cmd.color("dummy_70", "test2 and state 71", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 71")

cmd.set_color("dummy_71", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_71", "test2 and state 72", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 72")

cmd.set_color("dummy_72", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_72", "test2 and state 73", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 73")

cmd.set_color("dummy_73", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_73", "test2 and state 74", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 74")

cmd.set_color("dummy_74", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_74", "test2 and state 75", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 75")

cmd.set_color("dummy_75", (0.4971933871587852, 0.7122645136485968, 0.8380622837370243))
cmd.color("dummy_75", "test2 and state 76", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 76")

cmd.set_color("dummy_76", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_76", "test2 and state 77", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 77")

cmd.set_color("dummy_77", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_77", "test2 and state 78", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 78")

cmd.set_color("dummy_78", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_78", "test2 and state 79", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 79")

cmd.set_color("dummy_79", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_79", "test2 and state 80", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 80")

cmd.set_color("dummy_80", (0.7316416762783547, 0.08119953863898521, 0.15071126489811612))
cmd.color("dummy_80", "test2 and state 81", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 81")

cmd.set_color("dummy_81", (0.9951557093425606, 0.8322952710495963, 0.5213379469434832))
cmd.color("dummy_81", "test2 and state 82", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 82")

cmd.set_color("dummy_82", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_82", "test2 and state 83", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 83")

cmd.set_color("dummy_83", (0.993925413302576, 0.7707804690503652, 0.463514033064206))
cmd.color("dummy_83", "test2 and state 84", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 84")

cmd.set_color("dummy_84", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_84", "test2 and state 85", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 85")

cmd.set_color("dummy_85", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_85", "test2 and state 86", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 86")

cmd.set_color("dummy_86", (0.9078046905036524, 0.32425990003844674, 0.21537870049980778))
cmd.color("dummy_86", "test2 and state 87", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 87")

cmd.set_color("dummy_87", (0.2290657439446367, 0.3280276816608997, 0.641522491349481))
cmd.color("dummy_87", "test2 and state 88", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 88")

cmd.set_color("dummy_88", (0.368166089965398, 0.5750865051903116, 0.766089965397924))
cmd.color("dummy_88", "test2 and state 89", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 89")

cmd.set_color("dummy_89", (0.8632064590542099, 0.23044982698961938, 0.17231833910034605))
cmd.color("dummy_89", "test2 and state 90", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 90")

cmd.set_color("dummy_90", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_90", "test2 and state 91", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 91")

cmd.set_color("dummy_91", (0.9845444059976932, 0.6273740868896578, 0.3550173010380623))
cmd.color("dummy_91", "test2 and state 92", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 92")

cmd.set_color("dummy_92", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_92", "test2 and state 93", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 93")

cmd.set_color("dummy_93", (0.3898500576701269, 0.6009227220299886, 0.7794694348327567))
cmd.color("dummy_93", "test2 and state 94", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 94")

cmd.set_color("dummy_94", (0.9997693194925029, 0.9928489042675894, 0.7381776239907728))
cmd.color("dummy_94", "test2 and state 95", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 95")

cmd.set_color("dummy_95", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_95", "test2 and state 96", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 96")

cmd.set_color("dummy_96", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_96", "test2 and state 97", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 97")

cmd.set_color("dummy_97", (0.25059592464436753, 0.395847750865052, 0.6748942714340639))
cmd.color("dummy_97", "test2 and state 98", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 98")

cmd.set_color("dummy_98", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_98", "test2 and state 99", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 99")

cmd.set_color("dummy_99", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_99", "test2 and state 100", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 100")

cmd.set_color("dummy_100", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_100", "test2 and state 101", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 101")

cmd.set_color("dummy_101", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_101", "test2 and state 102", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 102")

cmd.set_color("dummy_102", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_102", "test2 and state 103", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 103")

cmd.set_color("dummy_103", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_103", "test2 and state 104", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 104")

cmd.set_color("dummy_104", (0.3898500576701269, 0.6009227220299886, 0.7794694348327567))
cmd.color("dummy_104", "test2 and state 105", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 105")

cmd.set_color("dummy_105", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_105", "test2 and state 106", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 106")

cmd.set_color("dummy_106", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_106", "test2 and state 107", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 107")

cmd.set_color("dummy_107", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_107", "test2 and state 108", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 108")

cmd.set_color("dummy_108", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_108", "test2 and state 109", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 109")

cmd.set_color("dummy_109", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_109", "test2 and state 110", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 110")

cmd.set_color("dummy_110", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_110", "test2 and state 111", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 111")

cmd.set_color("dummy_111", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_111", "test2 and state 112", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 112")

cmd.set_color("dummy_112", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_112", "test2 and state 113", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 113")

cmd.set_color("dummy_113", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_113", "test2 and state 114", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 114")

cmd.set_color("dummy_114", (0.9965397923875433, 0.8927335640138409, 0.5863898500576701))
cmd.color("dummy_114", "test2 and state 115", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 115")

cmd.set_color("dummy_115", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_115", "test2 and state 116", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 116")

cmd.set_color("dummy_116", (0.9933102652825836, 0.7400230680507497, 0.43460207612456747))
cmd.color("dummy_116", "test2 and state 117", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 117")

cmd.set_color("dummy_117", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_117", "test2 and state 118", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 118")

cmd.set_color("dummy_118", (0.5648596693579393, 0.7663975394079201, 0.8675893886966551))
cmd.color("dummy_118", "test2 and state 119", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 119")

cmd.set_color("dummy_119", (0.9996155324875048, 0.988081507112649, 0.7309496347558632))
cmd.color("dummy_119", "test2 and state 120", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 120")

cmd.set_color("dummy_120", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_120", "test2 and state 121", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 121")

cmd.set_color("dummy_121", (0.9665513264129182, 0.4974240676662822, 0.295040369088812))
cmd.color("dummy_121", "test2 and state 122", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 122")

cmd.set_color("dummy_122", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_122", "test2 and state 123", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 123")

cmd.set_color("dummy_123", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_123", "test2 and state 124", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 124")

cmd.set_color("dummy_124", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_124", "test2 and state 125", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 125")

cmd.set_color("dummy_125", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_125", "test2 and state 126", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 126")

cmd.set_color("dummy_126", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_126", "test2 and state 127", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 127")

cmd.set_color("dummy_127", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_127", "test2 and state 128", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 128")

cmd.set_color("dummy_128", (0.9993079584775086, 0.9785467128027683, 0.716493656286044))
cmd.color("dummy_128", "test2 and state 129", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 129")

cmd.set_color("dummy_129", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_129", "test2 and state 130", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 130")

cmd.set_color("dummy_130", (0.5564013840830451, 0.7596309111880047, 0.8638985005767013))
cmd.color("dummy_130", "test2 and state 131", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 131")

cmd.set_color("dummy_131", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_131", "test2 and state 132", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 132")

cmd.set_color("dummy_132", (0.5310265282583623, 0.7393310265282584, 0.8528258362168397))
cmd.color("dummy_132", "test2 and state 133", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 133")

cmd.set_color("dummy_133", (0.9833141099577086, 0.9935409457900808, 0.7797001153402537))
cmd.color("dummy_133", "test2 and state 134", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 134")

cmd.set_color("dummy_134", (0.8975009611687813, 0.960322952710496, 0.9374855824682813))
cmd.color("dummy_134", "test2 and state 135", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 135")

cmd.set_color("dummy_135", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_135", "test2 and state 136", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 136")

cmd.set_color("dummy_136", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_136", "test2 and state 137", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 137")

cmd.set_color("dummy_137", (0.47181853133410234, 0.6919646289888505, 0.8269896193771626))
cmd.color("dummy_137", "test2 and state 138", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 138")

cmd.set_color("dummy_138", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_138", "test2 and state 139", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 139")

cmd.set_color("dummy_139", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_139", "test2 and state 140", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 140")

cmd.set_color("dummy_140", (0.9928489042675894, 0.9972318339100346, 0.7621683967704729))
cmd.color("dummy_140", "test2 and state 141", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 141")

cmd.set_color("dummy_141", (0.7547097270280662, 0.10334486735870818, 0.15117262591311034))
cmd.color("dummy_141", "test2 and state 142", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 142")

cmd.set_color("dummy_142", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_142", "test2 and state 143", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 143")

cmd.set_color("dummy_143", (0.5141099577085737, 0.7257977700884276, 0.845444059976932))
cmd.color("dummy_143", "test2 and state 144", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 144")

cmd.set_color("dummy_144", (0.7623990772779701, 0.11072664359861592, 0.15132641291810842))
cmd.color("dummy_144", "test2 and state 145", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 145")

cmd.set_color("dummy_145", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_145", "test2 and state 146", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 146")

cmd.set_color("dummy_146", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_146", "test2 and state 147", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 147")

cmd.set_color("dummy_147", (0.9900807381776241, 0.6673587081891583, 0.3734717416378317))
cmd.color("dummy_147", "test2 and state 148", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 148")

cmd.set_color("dummy_148", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_148", "test2 and state 149", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 149")

cmd.set_color("dummy_149", (0.9637831603229527, 0.47743175701653207, 0.28581314878892733))
cmd.color("dummy_149", "test2 and state 150", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 150")

cmd.set_color("dummy_150", (0.9976163014225298, 0.9990772779700116, 0.7534025374855824))
cmd.color("dummy_150", "test2 and state 151", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 151")

cmd.set_color("dummy_151", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_151", "test2 and state 152", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 152")

cmd.set_color("dummy_152", (0.9996155324875048, 0.988081507112649, 0.7309496347558632))
cmd.color("dummy_152", "test2 and state 153", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 153")

cmd.set_color("dummy_153", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_153", "test2 and state 154", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 154")

cmd.set_color("dummy_154", (0.3537101114955786, 0.5578623606305267, 0.7571703191080355))
cmd.color("dummy_154", "test2 and state 155", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 155")

cmd.set_color("dummy_155", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_155", "test2 and state 156", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 156")

cmd.set_color("dummy_156", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_156", "test2 and state 157", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 157")

cmd.set_color("dummy_157", (0.6240676662821992, 0.8137639369473281, 0.8934256055363322))
cmd.color("dummy_157", "test2 and state 158", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 158")

cmd.set_color("dummy_158", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_158", "test2 and state 159", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 159")

cmd.set_color("dummy_159", (0.9959246443675509, 0.8707420222991157, 0.5574778931180315))
cmd.color("dummy_159", "test2 and state 160", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 160")

cmd.set_color("dummy_160", (0.9737793156478277, 0.9898500576701268, 0.7972318339100347))
cmd.color("dummy_160", "test2 and state 161", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 161")

cmd.set_color("dummy_161", (0.8702806612841217, 0.9489427143406383, 0.9702422145328721))
cmd.color("dummy_161", "test2 and state 162", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 162")

cmd.set_color("dummy_162", (0.5310265282583623, 0.7393310265282584, 0.8528258362168397))
cmd.color("dummy_162", "test2 and state 163", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 163")

cmd.set_color("dummy_163", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_163", "test2 and state 164", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 164")

cmd.set_color("dummy_164", (0.9968473663975395, 0.9022683583237218, 0.6008458285274896))
cmd.color("dummy_164", "test2 and state 165", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 165")

cmd.set_color("dummy_165", (0.8702806612841217, 0.9489427143406383, 0.9702422145328721))
cmd.color("dummy_165", "test2 and state 166", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 166")

cmd.set_color("dummy_166", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_166", "test2 and state 167", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 167")

cmd.set_color("dummy_167", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_167", "test2 and state 168", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 168")

cmd.set_color("dummy_168", (0.45490196078431383, 0.6784313725490198, 0.8196078431372549))
cmd.color("dummy_168", "test2 and state 169", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 169")

cmd.set_color("dummy_169", (0.9943867743175702, 0.7938485198000771, 0.4851980007689352))
cmd.color("dummy_169", "test2 and state 170", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 170")

cmd.set_color("dummy_170", (0.9833141099577086, 0.9935409457900808, 0.7797001153402537))
cmd.color("dummy_170", "test2 and state 171", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 171")

cmd.set_color("dummy_171", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_171", "test2 and state 172", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 172")

cmd.set_color("dummy_172", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_172", "test2 and state 173", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 173")

cmd.set_color("dummy_173", (0.4887351018838909, 0.7054978854286814, 0.8343713956170704))
cmd.color("dummy_173", "test2 and state 174", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 174")

cmd.set_color("dummy_174", (0.7724721261053442, 0.900961168781238, 0.942560553633218))
cmd.color("dummy_174", "test2 and state 175", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 175")

cmd.set_color("dummy_175", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_175", "test2 and state 176", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 176")

cmd.set_color("dummy_176", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_176", "test2 and state 177", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 177")

cmd.set_color("dummy_177", (0.9707035755478662, 0.5274125336409073, 0.308881199538639))
cmd.color("dummy_177", "test2 and state 178", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 178")

cmd.set_color("dummy_178", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_178", "test2 and state 179", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 179")

cmd.set_color("dummy_179", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_179", "test2 and state 180", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 180")

cmd.set_color("dummy_180", (0.952402921953095, 0.4180699730872741, 0.2584390618992695))
cmd.color("dummy_180", "test2 and state 181", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 181")

cmd.set_color("dummy_181", (0.5479430988081508, 0.7528642829680893, 0.8602076124567475))
cmd.color("dummy_181", "test2 and state 182", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 182")

cmd.set_color("dummy_182", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_182", "test2 and state 183", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 183")

cmd.set_color("dummy_183", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_183", "test2 and state 184", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 184")

cmd.set_color("dummy_184", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_184", "test2 and state 185", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 185")

cmd.set_color("dummy_185", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_185", "test2 and state 186", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 186")

cmd.set_color("dummy_186", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_186", "test2 and state 187", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 187")

cmd.set_color("dummy_187", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_187", "test2 and state 188", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 188")

cmd.set_color("dummy_188", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_188", "test2 and state 189", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 189")

cmd.set_color("dummy_189", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_189", "test2 and state 190", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 190")

cmd.set_color("dummy_190", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_190", "test2 and state 191", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 191")

cmd.set_color("dummy_191", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_191", "test2 and state 192", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 192")

cmd.set_color("dummy_192", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_192", "test2 and state 193", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 193")

cmd.set_color("dummy_193", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_193", "test2 and state 194", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 194")

cmd.set_color("dummy_194", (0.28143021914648214, 0.4717416378316033, 0.7125720876585929))
cmd.color("dummy_194", "test2 and state 195", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 195")

cmd.set_color("dummy_195", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_195", "test2 and state 196", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 196")

cmd.set_color("dummy_196", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_196", "test2 and state 197", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 197")

cmd.set_color("dummy_197", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_197", "test2 and state 198", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 198")

cmd.set_color("dummy_198", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_198", "test2 and state 199", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 199")

cmd.set_color("dummy_199", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_199", "test2 and state 200", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 200")

cmd.set_color("dummy_200", (0.4404459823144945, 0.661207227989235, 0.8106881968473664))
cmd.color("dummy_200", "test2 and state 201", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 201")

cmd.set_color("dummy_201", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_201", "test2 and state 202", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 202")

cmd.set_color("dummy_202", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_202", "test2 and state 203", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 203")

cmd.set_color("dummy_203", (0.9973087274125336, 0.9165705497885428, 0.6225297962322184))
cmd.color("dummy_203", "test2 and state 204", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 204")

cmd.set_color("dummy_204", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_204", "test2 and state 205", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 205")

cmd.set_color("dummy_205", (0.9946943483275663, 0.8092272202998847, 0.4996539792387543))
cmd.color("dummy_205", "test2 and state 206", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 206")

cmd.set_color("dummy_206", (0.9996155324875048, 0.988081507112649, 0.7309496347558632))
cmd.color("dummy_206", "test2 and state 207", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 207")

cmd.set_color("dummy_207", (0.5056516724336794, 0.7190311418685121, 0.8417531718569781))
cmd.color("dummy_207", "test2 and state 208", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 208")

cmd.set_color("dummy_208", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_208", "test2 and state 209", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 209")

cmd.set_color("dummy_209", (0.9070357554786621, 0.9640138408304498, 0.9199538638985004))
cmd.color("dummy_209", "test2 and state 210", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 210")

cmd.set_color("dummy_210", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_210", "test2 and state 211", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 211")

cmd.set_color("dummy_211", (0.4404459823144945, 0.661207227989235, 0.8106881968473664))
cmd.color("dummy_211", "test2 and state 212", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 212")

cmd.set_color("dummy_212", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_212", "test2 and state 213", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 213")

cmd.set_color("dummy_213", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_213", "test2 and state 214", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 214")

cmd.set_color("dummy_214", (0.9707035755478662, 0.5274125336409073, 0.308881199538639))
cmd.color("dummy_214", "test2 and state 215", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 215")

cmd.set_color("dummy_215", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_215", "test2 and state 216", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 216")

cmd.set_color("dummy_216", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_216", "test2 and state 217", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 217")

cmd.set_color("dummy_217", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_217", "test2 and state 218", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 218")

cmd.set_color("dummy_218", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_218", "test2 and state 219", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 219")

cmd.set_color("dummy_219", (0.41153402537485584, 0.6267589388696656, 0.7928489042675894))
cmd.color("dummy_219", "test2 and state 220", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 220")

cmd.set_color("dummy_220", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_220", "test2 and state 221", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 221")

cmd.set_color("dummy_221", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_221", "test2 and state 222", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 222")

cmd.set_color("dummy_222", (0.9936178392925799, 0.7554017685505575, 0.44905805459438675))
cmd.color("dummy_222", "test2 and state 223", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 223")

cmd.set_color("dummy_223", (0.5986928104575164, 0.7934640522875818, 0.8823529411764706))
cmd.color("dummy_223", "test2 and state 224", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 224")

cmd.set_color("dummy_224", (0.9954632833525567, 0.847673971549404, 0.5357939254133026))
cmd.color("dummy_224", "test2 and state 225", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 225")

cmd.set_color("dummy_225", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_225", "test2 and state 226", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 226")

cmd.set_color("dummy_226", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_226", "test2 and state 227", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 227")

cmd.set_color("dummy_227", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_227", "test2 and state 228", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 228")

cmd.set_color("dummy_228", (0.4802768166089966, 0.6987312572087659, 0.8306805074971165))
cmd.color("dummy_228", "test2 and state 229", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 229")

cmd.set_color("dummy_229", (0.9985390234525182, 0.9547097270280661, 0.6803537101114956))
cmd.color("dummy_229", "test2 and state 230", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 230")

cmd.set_color("dummy_230", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_230", "test2 and state 231", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 231")

cmd.set_color("dummy_231", (0.9997693194925029, 0.9928489042675894, 0.7381776239907728))
cmd.color("dummy_231", "test2 and state 232", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 232")

cmd.set_color("dummy_232", (0.9991541714725106, 0.9737793156478277, 0.7092656670511341))
cmd.color("dummy_232", "test2 and state 233", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 233")

cmd.set_color("dummy_233", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_233", "test2 and state 234", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 234")

cmd.set_color("dummy_234", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_234", "test2 and state 235", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 235")

cmd.set_color("dummy_235", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_235", "test2 and state 236", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 236")

cmd.set_color("dummy_236", (0.5310265282583623, 0.7393310265282584, 0.8528258362168397))
cmd.color("dummy_236", "test2 and state 237", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 237")

cmd.set_color("dummy_237", (0.9390234525182622, 0.3899269511726259, 0.245520953479431))
cmd.color("dummy_237", "test2 and state 238", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 238")

cmd.set_color("dummy_238", (0.976239907727797, 0.5673971549404075, 0.3273356401384083))
cmd.color("dummy_238", "test2 and state 239", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 239")

cmd.set_color("dummy_239", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_239", "test2 and state 240", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 240")

cmd.set_color("dummy_240", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_240", "test2 and state 241", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 241")

cmd.set_color("dummy_241", (0.9637831603229527, 0.47743175701653207, 0.28581314878892733))
cmd.color("dummy_241", "test2 and state 242", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 242")

cmd.set_color("dummy_242", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_242", "test2 and state 243", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 243")

cmd.set_color("dummy_243", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_243", "test2 and state 244", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 244")

cmd.set_color("dummy_244", (0.39707804690503656, 0.6095347943098809, 0.783929257977701))
cmd.color("dummy_244", "test2 and state 245", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 245")

cmd.set_color("dummy_245", (0.29588619761630147, 0.488965782391388, 0.7214917339484814))
cmd.color("dummy_245", "test2 and state 246", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 246")

cmd.set_color("dummy_246", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_246", "test2 and state 247", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 247")

cmd.set_color("dummy_247", (0.9976163014225298, 0.9261053440984237, 0.6369857747020377))
cmd.color("dummy_247", "test2 and state 248", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 248")

cmd.set_color("dummy_248", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_248", "test2 and state 249", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 249")

cmd.set_color("dummy_249", (0.996078431372549, 0.8784313725490196, 0.5647058823529412))
cmd.color("dummy_249", "test2 and state 250", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 250")

cmd.set_color("dummy_250", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_250", "test2 and state 251", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 251")

cmd.set_color("dummy_251", (0.8676662821991542, 0.2398308342945021, 0.1766243752402922))
cmd.color("dummy_251", "test2 and state 252", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 252")

cmd.set_color("dummy_252", (0.5225682429834679, 0.7325643983083431, 0.8491349480968858))
cmd.color("dummy_252", "test2 and state 253", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 253")

cmd.set_color("dummy_253", (0.994079200307574, 0.7784698193002692, 0.4707420222991157))
cmd.color("dummy_253", "test2 and state 254", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 254")

cmd.set_color("dummy_254", (0.9345636293733179, 0.38054594386774315, 0.24121491733948483))
cmd.color("dummy_254", "test2 and state 255", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 255")

cmd.set_color("dummy_255", (0.9990003844675125, 0.9690119184928874, 0.7020376778162245))
cmd.color("dummy_255", "test2 and state 256", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 256")

cmd.set_color("dummy_256", (0.9965397923875433, 0.8927335640138409, 0.5863898500576701))
cmd.color("dummy_256", "test2 and state 257", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 257")

cmd.set_color("dummy_257", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_257", "test2 and state 258", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 258")

cmd.set_color("dummy_258", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_258", "test2 and state 259", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 259")

cmd.set_color("dummy_259", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_259", "test2 and state 260", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 260")

cmd.set_color("dummy_260", (0.40430603613994626, 0.6181468665897734, 0.7883890811226452))
cmd.color("dummy_260", "test2 and state 261", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 261")

cmd.set_color("dummy_261", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_261", "test2 and state 262", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 262")

cmd.set_color("dummy_262", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_262", "test2 and state 263", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 263")

cmd.set_color("dummy_263", (0.9991541714725106, 0.9737793156478277, 0.7092656670511341))
cmd.color("dummy_263", "test2 and state 264", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 264")

cmd.set_color("dummy_264", (0.4404459823144945, 0.661207227989235, 0.8106881968473664))
cmd.color("dummy_264", "test2 and state 265", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 265")

cmd.set_color("dummy_265", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_265", "test2 and state 266", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 266")

cmd.set_color("dummy_266", (0.9301038062283737, 0.3711649365628604, 0.23690888119953865))
cmd.color("dummy_266", "test2 and state 267", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 267")

cmd.set_color("dummy_267", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_267", "test2 and state 268", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 268")

cmd.set_color("dummy_268", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_268", "test2 and state 269", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 269")

cmd.set_color("dummy_269", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_269", "test2 and state 270", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 270")

cmd.set_color("dummy_270", (0.9213379469434834, 0.9695501730103806, 0.8936562860438291))
cmd.color("dummy_270", "test2 and state 271", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 271")

cmd.set_color("dummy_271", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_271", "test2 and state 272", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 272")

cmd.set_color("dummy_272", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_272", "test2 and state 273", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 273")

cmd.set_color("dummy_273", (0.34648212226066905, 0.5492502883506345, 0.7527104959630911))
cmd.color("dummy_273", "test2 and state 274", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 274")

cmd.set_color("dummy_274", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_274", "test2 and state 275", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 275")

cmd.set_color("dummy_275", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_275", "test2 and state 276", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 276")

cmd.set_color("dummy_276", (0.4187620146097656, 0.635371011149558, 0.7973087274125337))
cmd.color("dummy_276", "test2 and state 277", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 277")

cmd.set_color("dummy_277", (0.9986928104575163, 0.9594771241830066, 0.6875816993464052))
cmd.color("dummy_277", "test2 and state 278", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 278")

cmd.set_color("dummy_278", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_278", "test2 and state 279", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 279")

cmd.set_color("dummy_279", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_279", "test2 and state 280", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 280")

cmd.set_color("dummy_280", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_280", "test2 and state 281", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 281")

cmd.set_color("dummy_281", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_281", "test2 and state 282", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 282")

cmd.set_color("dummy_282", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_282", "test2 and state 283", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 283")

cmd.set_color("dummy_283", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_283", "test2 and state 284", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 284")

cmd.set_color("dummy_284", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_284", "test2 and state 285", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 285")

cmd.set_color("dummy_285", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_285", "test2 and state 286", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 286")

cmd.set_color("dummy_286", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_286", "test2 and state 287", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 287")

cmd.set_color("dummy_287", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_287", "test2 and state 288", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 288")

cmd.set_color("dummy_288", (0.9790080738177624, 0.5873894655901576, 0.336562860438293))
cmd.color("dummy_288", "test2 and state 289", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 289")

cmd.set_color("dummy_289", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_289", "test2 and state 290", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 290")

cmd.set_color("dummy_290", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_290", "test2 and state 291", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 291")

cmd.set_color("dummy_291", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_291", "test2 and state 292", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 292")

cmd.set_color("dummy_292", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_292", "test2 and state 293", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 293")

cmd.set_color("dummy_293", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_293", "test2 and state 294", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 294")

cmd.set_color("dummy_294", (0.9499423298731258, 0.9806228373702423, 0.8410611303344866))
cmd.color("dummy_294", "test2 and state 295", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 295")

cmd.set_color("dummy_295", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_295", "test2 and state 296", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 296")

cmd.set_color("dummy_296", (0.9547097270280662, 0.9824682814302191, 0.8322952710495962))
cmd.color("dummy_296", "test2 and state 297", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 297")

cmd.set_color("dummy_297", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_297", "test2 and state 298", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 298")

cmd.set_color("dummy_298", (0.9833141099577086, 0.9935409457900808, 0.7797001153402537))
cmd.color("dummy_298", "test2 and state 299", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 299")

cmd.set_color("dummy_299", (0.9988465974625145, 0.9642445213379469, 0.6948096885813149))
cmd.color("dummy_299", "test2 and state 300", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 300")

cmd.set_color("dummy_300", (0.9547097270280662, 0.9824682814302191, 0.8322952710495962))
cmd.color("dummy_300", "test2 and state 301", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 301")

cmd.set_color("dummy_301", (0.8831987697039602, 0.9547866205305652, 0.9637831603229524))
cmd.color("dummy_301", "test2 and state 302", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 302")

cmd.set_color("dummy_302", (0.6240676662821992, 0.8137639369473281, 0.8934256055363322))
cmd.color("dummy_302", "test2 and state 303", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 303")

cmd.set_color("dummy_303", (0.8539792387543255, 0.9409457900807382, 0.9656286043829296))
cmd.color("dummy_303", "test2 and state 304", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 304")

cmd.set_color("dummy_304", (0.9070357554786621, 0.9640138408304498, 0.9199538638985004))
cmd.color("dummy_304", "test2 and state 305", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 305")

cmd.set_color("dummy_305", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_305", "test2 and state 306", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 306")

cmd.set_color("dummy_306", (0.9434832756632064, 0.39930795847750866, 0.24982698961937716))
cmd.color("dummy_306", "test2 and state 307", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 307")

cmd.set_color("dummy_307", (0.9785467128027682, 0.9916955017301038, 0.7884659746251441))
cmd.color("dummy_307", "test2 and state 308", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 308")

cmd.set_color("dummy_308", (0.9479430988081508, 0.4086889657823914, 0.25413302575932334))
cmd.color("dummy_308", "test2 and state 309", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 309")

cmd.set_color("dummy_309", (0.9996155324875048, 0.988081507112649, 0.7309496347558632))
cmd.color("dummy_309", "test2 and state 310", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 310")

cmd.set_color("dummy_310", (0.5986928104575164, 0.7934640522875818, 0.8823529411764706))
cmd.color("dummy_310", "test2 and state 311", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 311")

cmd.set_color("dummy_311", (0.9965397923875433, 0.8927335640138409, 0.5863898500576701))
cmd.color("dummy_311", "test2 and state 312", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 312")

cmd.set_color("dummy_312", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_312", "test2 and state 313", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 313")

cmd.set_color("dummy_313", (0.8899653979238754, 0.28673587081891583, 0.19815455594002307))
cmd.color("dummy_313", "test2 and state 314", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 314")

cmd.set_color("dummy_314", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_314", "test2 and state 315", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 315")

cmd.set_color("dummy_315", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_315", "test2 and state 316", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 316")

cmd.set_color("dummy_316", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_316", "test2 and state 317", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 317")

cmd.set_color("dummy_317", (0.9070357554786621, 0.9640138408304498, 0.9199538638985004))
cmd.color("dummy_317", "test2 and state 318", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 318")

cmd.set_color("dummy_318", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_318", "test2 and state 319", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 319")

cmd.set_color("dummy_319", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_319", "test2 and state 320", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 320")

cmd.set_color("dummy_320", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_320", "test2 and state 321", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 321")

cmd.set_color("dummy_321", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_321", "test2 and state 322", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 322")

cmd.set_color("dummy_322", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_322", "test2 and state 323", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 323")

cmd.set_color("dummy_323", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_323", "test2 and state 324", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 324")

cmd.set_color("dummy_324", (0.9845444059976932, 0.6273740868896578, 0.3550173010380623))
cmd.color("dummy_324", "test2 and state 325", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 325")

cmd.set_color("dummy_325", (0.9679354094579009, 0.5074202229911575, 0.2996539792387545))
cmd.color("dummy_325", "test2 and state 326", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 326")

cmd.set_color("dummy_326", (0.9859284890426759, 0.6373702422145329, 0.3596309111880046))
cmd.color("dummy_326", "test2 and state 327", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 327")

cmd.set_color("dummy_327", (0.9963860053825452, 0.8879661668589004, 0.5791618608227604))
cmd.color("dummy_327", "test2 and state 328", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 328")

cmd.set_color("dummy_328", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_328", "test2 and state 329", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 329")

cmd.set_color("dummy_329", (0.9070357554786621, 0.9640138408304498, 0.9199538638985004))
cmd.color("dummy_329", "test2 and state 330", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 330")

cmd.set_color("dummy_330", (0.9776239907727797, 0.5773933102652826, 0.33194925028835065))
cmd.color("dummy_330", "test2 and state 331", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 331")

cmd.set_color("dummy_331", (0.9451749327181853, 0.9787773933102653, 0.8498269896193771))
cmd.color("dummy_331", "test2 and state 332", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 332")

cmd.set_color("dummy_332", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_332", "test2 and state 333", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 333")

cmd.set_color("dummy_333", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_333", "test2 and state 334", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 334")

cmd.set_color("dummy_334", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_334", "test2 and state 335", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 335")

cmd.set_color("dummy_335", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_335", "test2 and state 336", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 336")

cmd.set_color("dummy_336", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_336", "test2 and state 337", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 337")

cmd.set_color("dummy_337", (0.38262206843521723, 0.5923106497500962, 0.7750096116878125))
cmd.color("dummy_337", "test2 and state 338", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 338")

cmd.set_color("dummy_338", (0.9937716262975779, 0.7630911188004613, 0.4562860438292964))
cmd.color("dummy_338", "test2 and state 339", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 339")

cmd.set_color("dummy_339", (0.9990003844675125, 0.9690119184928874, 0.7020376778162245))
cmd.color("dummy_339", "test2 and state 340", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 340")

cmd.set_color("dummy_340", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_340", "test2 and state 341", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 341")

cmd.set_color("dummy_341", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_341", "test2 and state 342", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 342")

cmd.set_color("dummy_342", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_342", "test2 and state 343", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 343")

cmd.set_color("dummy_343", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_343", "test2 and state 344", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 344")

cmd.set_color("dummy_344", (0.6828143021914649, 0.8569780853517878, 0.9171856978085352))
cmd.color("dummy_344", "test2 and state 345", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 345")

cmd.set_color("dummy_345", (0.38262206843521723, 0.5923106497500962, 0.7750096116878125))
cmd.color("dummy_345", "test2 and state 346", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 346")

cmd.set_color("dummy_346", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_346", "test2 and state 347", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 347")

cmd.set_color("dummy_347", (0.9499423298731258, 0.9806228373702423, 0.8410611303344866))
cmd.color("dummy_347", "test2 and state 348", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 348")

cmd.set_color("dummy_348", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_348", "test2 and state 349", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 349")

cmd.set_color("dummy_349", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_349", "test2 and state 350", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 350")

cmd.set_color("dummy_350", (0.9451749327181853, 0.9787773933102653, 0.8498269896193771))
cmd.color("dummy_350", "test2 and state 351", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 351")

cmd.set_color("dummy_351", (0.9976163014225298, 0.9990772779700116, 0.7534025374855824))
cmd.color("dummy_351", "test2 and state 352", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 352")

cmd.set_color("dummy_352", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_352", "test2 and state 353", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 353")

cmd.set_color("dummy_353", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_353", "test2 and state 354", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 354")

cmd.set_color("dummy_354", (0.9213379469434834, 0.9695501730103806, 0.8936562860438291))
cmd.color("dummy_354", "test2 and state 355", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 355")

cmd.set_color("dummy_355", (0.9930026912725874, 0.724644367550942, 0.4201460976547482))
cmd.color("dummy_355", "test2 and state 356", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 356")

cmd.set_color("dummy_356", (0.9996155324875048, 0.988081507112649, 0.7309496347558632))
cmd.color("dummy_356", "test2 and state 357", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 357")

cmd.set_color("dummy_357", (0.4187620146097656, 0.635371011149558, 0.7973087274125337))
cmd.color("dummy_357", "test2 and state 358", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 358")

cmd.set_color("dummy_358", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_358", "test2 and state 359", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 359")

cmd.set_color("dummy_359", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_359", "test2 and state 360", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 360")

cmd.set_color("dummy_360", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_360", "test2 and state 361", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 361")

cmd.set_color("dummy_361", (0.9977700884275279, 0.930872741253364, 0.6442137639369473))
cmd.color("dummy_361", "test2 and state 362", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 362")

cmd.set_color("dummy_362", (0.6828143021914649, 0.8569780853517878, 0.9171856978085352))
cmd.color("dummy_362", "test2 and state 363", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 363")

cmd.set_color("dummy_363", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_363", "test2 and state 364", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 364")

cmd.set_color("dummy_364", (0.9951557093425606, 0.8322952710495963, 0.5213379469434832))
cmd.color("dummy_364", "test2 and state 365", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 365")

cmd.set_color("dummy_365", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_365", "test2 and state 366", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 366")

cmd.set_color("dummy_366", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_366", "test2 and state 367", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 367")

cmd.set_color("dummy_367", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_367", "test2 and state 368", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 368")

cmd.set_color("dummy_368", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_368", "test2 and state 369", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 369")

cmd.set_color("dummy_369", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_369", "test2 and state 370", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 370")

cmd.set_color("dummy_370", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_370", "test2 and state 371", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 371")

cmd.set_color("dummy_371", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_371", "test2 and state 372", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 372")

cmd.set_color("dummy_372", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_372", "test2 and state 373", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 373")

cmd.set_color("dummy_373", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_373", "test2 and state 374", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 374")

cmd.set_color("dummy_374", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_374", "test2 and state 375", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 375")

cmd.set_color("dummy_375", (0.9914648212226067, 0.6773548635140331, 0.3780853517877739))
cmd.color("dummy_375", "test2 and state 376", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 376")

cmd.set_color("dummy_376", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_376", "test2 and state 377", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 377")

cmd.set_color("dummy_377", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_377", "test2 and state 378", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 378")

cmd.set_color("dummy_378", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_378", "test2 and state 379", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 379")

cmd.set_color("dummy_379", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_379", "test2 and state 380", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 380")

cmd.set_color("dummy_380", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_380", "test2 and state 381", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 381")

cmd.set_color("dummy_381", (0.4971933871587852, 0.7122645136485968, 0.8380622837370243))
cmd.color("dummy_381", "test2 and state 382", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 382")

cmd.set_color("dummy_382", (0.8975009611687813, 0.960322952710496, 0.9374855824682813))
cmd.color("dummy_382", "test2 and state 383", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 383")

cmd.set_color("dummy_383", (0.9928489042675894, 0.7169550173010381, 0.4129181084198385))
cmd.color("dummy_383", "test2 and state 384", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 384")

cmd.set_color("dummy_384", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_384", "test2 and state 385", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 385")

cmd.set_color("dummy_385", (0.5986928104575164, 0.7934640522875818, 0.8823529411764706))
cmd.color("dummy_385", "test2 and state 386", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 386")

cmd.set_color("dummy_386", (0.9499423298731258, 0.9806228373702423, 0.8410611303344866))
cmd.color("dummy_386", "test2 and state 387", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 387")

cmd.set_color("dummy_387", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_387", "test2 and state 388", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 388")

cmd.set_color("dummy_388", (0.9922337562475971, 0.6861976163014225, 0.3840061514801999))
cmd.color("dummy_388", "test2 and state 389", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 389")

cmd.set_color("dummy_389", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_389", "test2 and state 390", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 390")

cmd.set_color("dummy_390", (0.9925413302575933, 0.7015763168012303, 0.3984621299500192))
cmd.color("dummy_390", "test2 and state 391", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 391")

cmd.set_color("dummy_391", (0.4633602460592081, 0.6851980007689351, 0.8232987312572088))
cmd.color("dummy_391", "test2 and state 392", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 392")

cmd.set_color("dummy_392", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_392", "test2 and state 393", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 393")

cmd.set_color("dummy_393", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_393", "test2 and state 394", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 394")

cmd.set_color("dummy_394", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_394", "test2 and state 395", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 395")

cmd.set_color("dummy_395", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_395", "test2 and state 396", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 396")

cmd.set_color("dummy_396", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_396", "test2 and state 397", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 397")

cmd.set_color("dummy_397", (0.9211841599384853, 0.3524029219530952, 0.22829680891964643))
cmd.color("dummy_397", "test2 and state 398", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 398")

cmd.set_color("dummy_398", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_398", "test2 and state 399", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 399")

cmd.set_color("dummy_399", (0.9451749327181853, 0.9787773933102653, 0.8498269896193771))
cmd.color("dummy_399", "test2 and state 400", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 400")

cmd.set_color("dummy_400", (0.9993079584775086, 0.9785467128027683, 0.716493656286044))
cmd.color("dummy_400", "test2 and state 401", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 401")

cmd.set_color("dummy_401", (0.9931564782775856, 0.7323337178008458, 0.42737408688965783))
cmd.color("dummy_401", "test2 and state 402", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 402")

cmd.set_color("dummy_402", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_402", "test2 and state 403", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 403")

cmd.set_color("dummy_403", (0.37539407920030765, 0.5836985774702039, 0.7705497885428682))
cmd.color("dummy_403", "test2 and state 404", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 404")

cmd.set_color("dummy_404", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_404", "test2 and state 405", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 405")

cmd.set_color("dummy_405", (0.9985390234525182, 0.9547097270280661, 0.6803537101114956))
cmd.color("dummy_405", "test2 and state 406", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 406")

cmd.set_color("dummy_406", (0.9948481353325644, 0.8169165705497885, 0.506881968473664))
cmd.color("dummy_406", "test2 and state 407", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 407")

cmd.set_color("dummy_407", (0.9610149942329873, 0.457439446366782, 0.2765859284890427))
cmd.color("dummy_407", "test2 and state 408", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 408")

cmd.set_color("dummy_408", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_408", "test2 and state 409", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 409")

cmd.set_color("dummy_409", (0.47181853133410234, 0.6919646289888505, 0.8269896193771626))
cmd.color("dummy_409", "test2 and state 410", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 410")

cmd.set_color("dummy_410", (0.9078046905036524, 0.32425990003844674, 0.21537870049980778))
cmd.color("dummy_410", "test2 and state 411", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 411")

cmd.set_color("dummy_411", (0.5056516724336794, 0.7190311418685121, 0.8417531718569781))
cmd.color("dummy_411", "test2 and state 412", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 412")

cmd.set_color("dummy_412", (0.9991541714725106, 0.9737793156478277, 0.7092656670511341))
cmd.color("dummy_412", "test2 and state 413", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 413")

cmd.set_color("dummy_413", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_413", "test2 and state 414", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 414")

cmd.set_color("dummy_414", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_414", "test2 and state 415", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 415")

cmd.set_color("dummy_415", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_415", "test2 and state 416", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 416")

cmd.set_color("dummy_416", (0.9945405613225683, 0.8015378700499808, 0.4924259900038447))
cmd.color("dummy_416", "test2 and state 417", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 417")

cmd.set_color("dummy_417", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_417", "test2 and state 418", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 418")

cmd.set_color("dummy_418", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_418", "test2 and state 419", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 419")

cmd.set_color("dummy_419", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_419", "test2 and state 420", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 420")

cmd.set_color("dummy_420", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_420", "test2 and state 421", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 421")

cmd.set_color("dummy_421", (0.40430603613994626, 0.6181468665897734, 0.7883890811226452))
cmd.color("dummy_421", "test2 and state 422", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 422")

cmd.set_color("dummy_422", (0.8975009611687813, 0.960322952710496, 0.9374855824682813))
cmd.color("dummy_422", "test2 and state 423", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 423")

cmd.set_color("dummy_423", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_423", "test2 and state 424", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 424")

cmd.set_color("dummy_424", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_424", "test2 and state 425", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 425")

cmd.set_color("dummy_425", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_425", "test2 and state 426", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 426")

cmd.set_color("dummy_426", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_426", "test2 and state 427", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 427")

cmd.set_color("dummy_427", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_427", "test2 and state 428", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 428")

cmd.set_color("dummy_428", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_428", "test2 and state 429", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 429")

cmd.set_color("dummy_429", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_429", "test2 and state 430", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 430")

cmd.set_color("dummy_430", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_430", "test2 and state 431", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 431")

cmd.set_color("dummy_431", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_431", "test2 and state 432", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 432")

cmd.set_color("dummy_432", (0.4887351018838909, 0.7054978854286814, 0.8343713956170704))
cmd.color("dummy_432", "test2 and state 433", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 433")

cmd.set_color("dummy_433", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_433", "test2 and state 434", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 434")

cmd.set_color("dummy_434", (0.3320261437908497, 0.5320261437908498, 0.7437908496732026))
cmd.color("dummy_434", "test2 and state 435", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 435")

cmd.set_color("dummy_435", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_435", "test2 and state 436", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 436")

cmd.set_color("dummy_436", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_436", "test2 and state 437", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 437")

cmd.set_color("dummy_437", (0.9974625144175318, 0.9213379469434833, 0.629757785467128))
cmd.color("dummy_437", "test2 and state 438", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 438")

cmd.set_color("dummy_438", (0.4259900038446752, 0.6439830834294503, 0.801768550557478))
cmd.color("dummy_438", "test2 and state 439", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 439")

cmd.set_color("dummy_439", (0.8831987697039602, 0.9547866205305652, 0.9637831603229524))
cmd.color("dummy_439", "test2 and state 440", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 440")

cmd.set_color("dummy_440", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_440", "test2 and state 441", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 441")

cmd.set_color("dummy_441", (0.9993079584775086, 0.9785467128027683, 0.716493656286044))
cmd.color("dummy_441", "test2 and state 442", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 442")

cmd.set_color("dummy_442", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_442", "test2 and state 443", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 443")

cmd.set_color("dummy_443", (0.43321799307958486, 0.6525951557093427, 0.8062283737024222))
cmd.color("dummy_443", "test2 and state 444", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 444")

cmd.set_color("dummy_444", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_444", "test2 and state 445", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 445")

cmd.set_color("dummy_445", (0.2659746251441753, 0.4442906574394464, 0.6987312572087659))
cmd.color("dummy_445", "test2 and state 446", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 446")

cmd.set_color("dummy_446", (0.9988465974625145, 0.9642445213379469, 0.6948096885813149))
cmd.color("dummy_446", "test2 and state 447", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 447")

cmd.set_color("dummy_447", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_447", "test2 and state 448", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 448")

cmd.set_color("dummy_448", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_448", "test2 and state 449", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 449")

cmd.set_color("dummy_449", (0.7724721261053442, 0.900961168781238, 0.942560553633218))
cmd.color("dummy_449", "test2 and state 450", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 450")

cmd.set_color("dummy_450", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_450", "test2 and state 451", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 451")

cmd.set_color("dummy_451", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_451", "test2 and state 452", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 452")

cmd.set_color("dummy_452", (0.5648596693579393, 0.7663975394079201, 0.8675893886966551))
cmd.color("dummy_452", "test2 and state 453", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 453")

cmd.set_color("dummy_453", (0.958246828143022, 0.43744713571703187, 0.267358708189158))
cmd.color("dummy_453", "test2 and state 454", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 454")

cmd.set_color("dummy_454", (0.9965397923875433, 0.8927335640138409, 0.5863898500576701))
cmd.color("dummy_454", "test2 and state 455", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 455")

cmd.set_color("dummy_455", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_455", "test2 and state 456", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 456")

cmd.set_color("dummy_456", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_456", "test2 and state 457", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 457")

cmd.set_color("dummy_457", (0.9962322183775472, 0.88319876970396, 0.5719338715878508))
cmd.color("dummy_457", "test2 and state 458", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 458")

cmd.set_color("dummy_458", (0.5479430988081508, 0.7528642829680893, 0.8602076124567475))
cmd.color("dummy_458", "test2 and state 459", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 459")

cmd.set_color("dummy_459", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_459", "test2 and state 460", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 460")

cmd.set_color("dummy_460", (0.47181853133410234, 0.6919646289888505, 0.8269896193771626))
cmd.color("dummy_460", "test2 and state 461", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 461")

cmd.set_color("dummy_461", (0.9937716262975779, 0.7630911188004613, 0.4562860438292964))
cmd.color("dummy_461", "test2 and state 462", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 462")

cmd.set_color("dummy_462", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_462", "test2 and state 463", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 463")

cmd.set_color("dummy_463", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_463", "test2 and state 464", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 464")

cmd.set_color("dummy_464", (0.9993079584775086, 0.9785467128027683, 0.716493656286044))
cmd.color("dummy_464", "test2 and state 465", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 465")

cmd.set_color("dummy_465", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_465", "test2 and state 466", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 466")

cmd.set_color("dummy_466", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_466", "test2 and state 467", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 467")

cmd.set_color("dummy_467", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_467", "test2 and state 468", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 468")

cmd.set_color("dummy_468", (0.9985390234525182, 0.9547097270280661, 0.6803537101114956))
cmd.color("dummy_468", "test2 and state 469", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 469")

cmd.set_color("dummy_469", (0.9737793156478277, 0.9898500576701268, 0.7972318339100347))
cmd.color("dummy_469", "test2 and state 470", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 470")

cmd.set_color("dummy_470", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_470", "test2 and state 471", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 471")

cmd.set_color("dummy_471", (0.9953094963475586, 0.8399846212995001, 0.5285659361783929))
cmd.color("dummy_471", "test2 and state 472", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 472")

cmd.set_color("dummy_472", (0.4259900038446752, 0.6439830834294503, 0.801768550557478))
cmd.color("dummy_472", "test2 and state 473", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 473")

cmd.set_color("dummy_473", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_473", "test2 and state 474", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 474")

cmd.set_color("dummy_474", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_474", "test2 and state 475", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 475")

cmd.set_color("dummy_475", (0.9845444059976932, 0.6273740868896578, 0.3550173010380623))
cmd.color("dummy_475", "test2 and state 476", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 476")

cmd.set_color("dummy_476", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_476", "test2 and state 477", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 477")

cmd.set_color("dummy_477", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_477", "test2 and state 478", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 478")

cmd.set_color("dummy_478", (0.9974625144175318, 0.9213379469434833, 0.629757785467128))
cmd.color("dummy_478", "test2 and state 479", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 479")

cmd.set_color("dummy_479", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_479", "test2 and state 480", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 480")

cmd.set_color("dummy_480", (0.9928489042675894, 0.7169550173010381, 0.4129181084198385))
cmd.color("dummy_480", "test2 and state 481", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 481")

cmd.set_color("dummy_481", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_481", "test2 and state 482", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 482")

cmd.set_color("dummy_482", (0.9979238754325259, 0.9356401384083045, 0.6514417531718569))
cmd.color("dummy_482", "test2 and state 483", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 483")

cmd.set_color("dummy_483", (0.996078431372549, 0.8784313725490196, 0.5647058823529412))
cmd.color("dummy_483", "test2 and state 484", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 484")

cmd.set_color("dummy_484", (0.4971933871587852, 0.7122645136485968, 0.8380622837370243))
cmd.color("dummy_484", "test2 and state 485", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 485")

cmd.set_color("dummy_485", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_485", "test2 and state 486", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 486")

cmd.set_color("dummy_486", (0.996078431372549, 0.8784313725490196, 0.5647058823529412))
cmd.color("dummy_486", "test2 and state 487", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 487")

cmd.set_color("dummy_487", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_487", "test2 and state 488", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 488")

cmd.set_color("dummy_488", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_488", "test2 and state 489", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 489")

cmd.set_color("dummy_489", (0.8975009611687813, 0.960322952710496, 0.9374855824682813))
cmd.color("dummy_489", "test2 and state 490", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 490")

cmd.set_color("dummy_490", (0.9982314494425221, 0.9451749327181853, 0.6658977316416763))
cmd.color("dummy_490", "test2 and state 491", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 491")

cmd.set_color("dummy_491", (0.9976163014225298, 0.9990772779700116, 0.7534025374855824))
cmd.color("dummy_491", "test2 and state 492", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 492")

cmd.set_color("dummy_492", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_492", "test2 and state 493", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 493")

cmd.set_color("dummy_493", (0.32479815455594, 0.5234140715109573, 0.7393310265282584))
cmd.color("dummy_493", "test2 and state 494", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 494")

cmd.set_color("dummy_494", (0.5056516724336794, 0.7190311418685121, 0.8417531718569781))
cmd.color("dummy_494", "test2 and state 495", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 495")

cmd.set_color("dummy_495", (0.21676278354479048, 0.2892733564013841, 0.6224529027297194))
cmd.color("dummy_495", "test2 and state 496", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 496")

cmd.set_color("dummy_496", (0.9390234525182622, 0.3899269511726259, 0.245520953479431))
cmd.color("dummy_496", "test2 and state 497", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 497")

cmd.set_color("dummy_497", (0.9957708573625529, 0.8630526720492119, 0.5502499038831219))
cmd.color("dummy_497", "test2 and state 498", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 498")

cmd.set_color("dummy_498", (0.28143021914648214, 0.4717416378316033, 0.7125720876585929))
cmd.color("dummy_498", "test2 and state 499", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 499")

cmd.set_color("dummy_499", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_499", "test2 and state 500", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 500")

cmd.set_color("dummy_500", (0.6470588235294118, 0.0, 0.14901960784313725))
cmd.color("dummy_500", "test2 and state 501", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 501")

cmd.set_color("dummy_501", (0.9943867743175702, 0.7938485198000771, 0.4851980007689352))
cmd.color("dummy_501", "test2 and state 502", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 502")

cmd.set_color("dummy_502", (0.45490196078431383, 0.6784313725490198, 0.8196078431372549))
cmd.color("dummy_502", "test2 and state 503", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 503")

cmd.set_color("dummy_503", (0.9737793156478277, 0.9898500576701268, 0.7972318339100347))
cmd.color("dummy_503", "test2 and state 504", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 504")

cmd.set_color("dummy_504", (0.6240676662821992, 0.8137639369473281, 0.8934256055363322))
cmd.color("dummy_504", "test2 and state 505", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 505")

cmd.set_color("dummy_505", (0.9665513264129182, 0.4974240676662822, 0.295040369088812))
cmd.color("dummy_505", "test2 and state 506", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 506")

cmd.set_color("dummy_506", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_506", "test2 and state 507", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 507")

cmd.set_color("dummy_507", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_507", "test2 and state 508", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 508")

cmd.set_color("dummy_508", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_508", "test2 and state 509", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 509")

cmd.set_color("dummy_509", (0.9956170703575548, 0.855363321799308, 0.5430219146482123))
cmd.color("dummy_509", "test2 and state 510", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 510")

cmd.set_color("dummy_510", (0.5648596693579393, 0.7663975394079201, 0.8675893886966551))
cmd.color("dummy_510", "test2 and state 511", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 511")

cmd.set_color("dummy_511", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_511", "test2 and state 512", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 512")

cmd.set_color("dummy_512", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_512", "test2 and state 513", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 513")

cmd.set_color("dummy_513", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_513", "test2 and state 514", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 514")

cmd.set_color("dummy_514", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_514", "test2 and state 515", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 515")

cmd.set_color("dummy_515", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_515", "test2 and state 516", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 516")

cmd.set_color("dummy_516", (0.9776239907727797, 0.5773933102652826, 0.33194925028835065))
cmd.color("dummy_516", "test2 and state 517", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 517")

cmd.set_color("dummy_517", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_517", "test2 and state 518", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 518")

cmd.set_color("dummy_518", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_518", "test2 and state 519", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 519")

cmd.set_color("dummy_519", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_519", "test2 and state 520", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 520")

cmd.set_color("dummy_520", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_520", "test2 and state 521", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 521")

cmd.set_color("dummy_521", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_521", "test2 and state 522", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 522")

cmd.set_color("dummy_522", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_522", "test2 and state 523", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 523")

cmd.set_color("dummy_523", (0.5479430988081508, 0.7528642829680893, 0.8602076124567475))
cmd.color("dummy_523", "test2 and state 524", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 524")

cmd.set_color("dummy_524", (0.9785467128027682, 0.9916955017301038, 0.7884659746251441))
cmd.color("dummy_524", "test2 and state 525", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 525")

cmd.set_color("dummy_525", (0.9990003844675125, 0.9690119184928874, 0.7020376778162245))
cmd.color("dummy_525", "test2 and state 526", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 526")

cmd.set_color("dummy_526", (0.9963860053825452, 0.8879661668589004, 0.5791618608227604))
cmd.color("dummy_526", "test2 and state 527", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 527")

cmd.set_color("dummy_527", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_527", "test2 and state 528", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 528")

cmd.set_color("dummy_528", (0.9997693194925029, 0.9928489042675894, 0.7381776239907728))
cmd.color("dummy_528", "test2 and state 529", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 529")

cmd.set_color("dummy_529", (0.9974625144175318, 0.9213379469434833, 0.629757785467128))
cmd.color("dummy_529", "test2 and state 530", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 530")

cmd.set_color("dummy_530", (0.9999231064975009, 0.9976163014225298, 0.7454056132256824))
cmd.color("dummy_530", "test2 and state 531", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 531")

cmd.set_color("dummy_531", (0.8975009611687813, 0.960322952710496, 0.9374855824682813))
cmd.color("dummy_531", "test2 and state 532", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 532")

cmd.set_color("dummy_532", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_532", "test2 and state 533", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 533")

cmd.set_color("dummy_533", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_533", "test2 and state 534", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 534")

cmd.set_color("dummy_534", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_534", "test2 and state 535", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 535")

cmd.set_color("dummy_535", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_535", "test2 and state 536", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 536")

cmd.set_color("dummy_536", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_536", "test2 and state 537", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 537")

cmd.set_color("dummy_537", (0.9942329873125721, 0.7861591695501731, 0.47797001153402535))
cmd.color("dummy_537", "test2 and state 538", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 538")

cmd.set_color("dummy_538", (0.9213379469434834, 0.9695501730103806, 0.8936562860438291))
cmd.color("dummy_538", "test2 and state 539", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 539")

cmd.set_color("dummy_539", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_539", "test2 and state 540", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 540")

cmd.set_color("dummy_540", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_540", "test2 and state 541", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 541")

cmd.set_color("dummy_541", (0.9167243367935409, 0.3430219146482122, 0.22399077277970011))
cmd.color("dummy_541", "test2 and state 542", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 542")

cmd.set_color("dummy_542", (0.8975009611687813, 0.960322952710496, 0.9374855824682813))
cmd.color("dummy_542", "test2 and state 543", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 543")

cmd.set_color("dummy_543", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_543", "test2 and state 544", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 544")

cmd.set_color("dummy_544", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_544", "test2 and state 545", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 545")

cmd.set_color("dummy_545", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_545", "test2 and state 546", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 546")

cmd.set_color("dummy_546", (0.9966935793925413, 0.8975009611687812, 0.5936178392925797))
cmd.color("dummy_546", "test2 and state 547", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 547")

cmd.set_color("dummy_547", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_547", "test2 and state 548", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 548")

cmd.set_color("dummy_548", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_548", "test2 and state 549", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 549")

cmd.set_color("dummy_549", (0.9951557093425606, 0.8322952710495963, 0.5213379469434832))
cmd.color("dummy_549", "test2 and state 550", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 550")

cmd.set_color("dummy_550", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_550", "test2 and state 551", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 551")

cmd.set_color("dummy_551", (0.9078046905036524, 0.32425990003844674, 0.21537870049980778))
cmd.color("dummy_551", "test2 and state 552", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 552")

cmd.set_color("dummy_552", (0.34648212226066905, 0.5492502883506345, 0.7527104959630911))
cmd.color("dummy_552", "test2 and state 553", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 553")

cmd.set_color("dummy_553", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_553", "test2 and state 554", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 554")

cmd.set_color("dummy_554", (0.40430603613994626, 0.6181468665897734, 0.7883890811226452))
cmd.color("dummy_554", "test2 and state 555", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 555")

cmd.set_color("dummy_555", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_555", "test2 and state 556", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 556")

cmd.set_color("dummy_556", (0.3898500576701269, 0.6009227220299886, 0.7794694348327567))
cmd.color("dummy_556", "test2 and state 557", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 557")

cmd.set_color("dummy_557", (0.8721261053440984, 0.24921184159938484, 0.18093041138023838))
cmd.color("dummy_557", "test2 and state 558", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 558")

cmd.set_color("dummy_558", (0.3320261437908497, 0.5320261437908498, 0.7437908496732026))
cmd.color("dummy_558", "test2 and state 559", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 559")

cmd.set_color("dummy_559", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_559", "test2 and state 560", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 560")

cmd.set_color("dummy_560", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_560", "test2 and state 561", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 561")

cmd.set_color("dummy_561", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_561", "test2 and state 562", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 562")

cmd.set_color("dummy_562", (0.9974625144175318, 0.9213379469434833, 0.629757785467128))
cmd.color("dummy_562", "test2 and state 563", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 563")

cmd.set_color("dummy_563", (0.9974625144175318, 0.9213379469434833, 0.629757785467128))
cmd.color("dummy_563", "test2 and state 564", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 564")

cmd.set_color("dummy_564", (0.9928489042675894, 0.7169550173010381, 0.4129181084198385))
cmd.color("dummy_564", "test2 and state 565", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 565")

cmd.set_color("dummy_565", (0.8810457516339869, 0.2679738562091503, 0.18954248366013074))
cmd.color("dummy_565", "test2 and state 566", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 566")

cmd.set_color("dummy_566", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_566", "test2 and state 567", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 567")

cmd.set_color("dummy_567", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_567", "test2 and state 568", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 568")

cmd.set_color("dummy_568", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_568", "test2 and state 569", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 569")

cmd.set_color("dummy_569", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_569", "test2 and state 570", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 570")

cmd.set_color("dummy_570", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_570", "test2 and state 571", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 571")

cmd.set_color("dummy_571", (0.4404459823144945, 0.661207227989235, 0.8106881968473664))
cmd.color("dummy_571", "test2 and state 572", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 572")

cmd.set_color("dummy_572", (0.9943867743175702, 0.7938485198000771, 0.4851980007689352))
cmd.color("dummy_572", "test2 and state 573", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 573")

cmd.set_color("dummy_573", (0.9971549404075356, 0.9118031526336025, 0.6153018069973087))
cmd.color("dummy_573", "test2 and state 574", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 574")

cmd.set_color("dummy_574", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_574", "test2 and state 575", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 575")

cmd.set_color("dummy_575", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_575", "test2 and state 576", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 576")

cmd.set_color("dummy_576", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_576", "test2 and state 577", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 577")

cmd.set_color("dummy_577", (0.7724721261053442, 0.900961168781238, 0.942560553633218))
cmd.color("dummy_577", "test2 and state 578", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 578")

cmd.set_color("dummy_578", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_578", "test2 and state 579", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 579")

cmd.set_color("dummy_579", (0.5986928104575164, 0.7934640522875818, 0.8823529411764706))
cmd.color("dummy_579", "test2 and state 580", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 580")

cmd.set_color("dummy_580", (0.9547097270280662, 0.9824682814302191, 0.8322952710495962))
cmd.color("dummy_580", "test2 and state 581", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 581")

cmd.set_color("dummy_581", (0.9831603229527105, 0.6173779315647828, 0.35040369088811996))
cmd.color("dummy_581", "test2 and state 582", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 582")

cmd.set_color("dummy_582", (0.41153402537485584, 0.6267589388696656, 0.7928489042675894))
cmd.color("dummy_582", "test2 and state 583", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 583")

cmd.set_color("dummy_583", (0.9963860053825452, 0.8879661668589004, 0.5791618608227604))
cmd.color("dummy_583", "test2 and state 584", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 584")

cmd.set_color("dummy_584", (0.9078046905036524, 0.32425990003844674, 0.21537870049980778))
cmd.color("dummy_584", "test2 and state 585", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 585")

cmd.set_color("dummy_585", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_585", "test2 and state 586", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 586")

cmd.set_color("dummy_586", (0.9996155324875048, 0.988081507112649, 0.7309496347558632))
cmd.color("dummy_586", "test2 and state 587", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 587")

cmd.set_color("dummy_587", (0.5225682429834679, 0.7325643983083431, 0.8491349480968858))
cmd.color("dummy_587", "test2 and state 588", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 588")

cmd.set_color("dummy_588", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_588", "test2 and state 589", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 589")

cmd.set_color("dummy_589", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_589", "test2 and state 590", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 590")

cmd.set_color("dummy_590", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_590", "test2 and state 591", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 591")

cmd.set_color("dummy_591", (0.5056516724336794, 0.7190311418685121, 0.8417531718569781))
cmd.color("dummy_591", "test2 and state 592", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 592")

cmd.set_color("dummy_592", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_592", "test2 and state 593", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 593")

cmd.set_color("dummy_593", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_593", "test2 and state 594", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 594")

cmd.set_color("dummy_594", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_594", "test2 and state 595", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 595")

cmd.set_color("dummy_595", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_595", "test2 and state 596", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 596")

cmd.set_color("dummy_596", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_596", "test2 and state 597", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 597")

cmd.set_color("dummy_597", (0.4259900038446752, 0.6439830834294503, 0.801768550557478))
cmd.color("dummy_597", "test2 and state 598", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 598")

cmd.set_color("dummy_598", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_598", "test2 and state 599", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 599")

cmd.set_color("dummy_599", (0.5056516724336794, 0.7190311418685121, 0.8417531718569781))
cmd.color("dummy_599", "test2 and state 600", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 600")

cmd.set_color("dummy_600", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_600", "test2 and state 601", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 601")

cmd.set_color("dummy_601", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_601", "test2 and state 602", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 602")

cmd.set_color("dummy_602", (0.40430603613994626, 0.6181468665897734, 0.7883890811226452))
cmd.color("dummy_602", "test2 and state 603", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 603")

cmd.set_color("dummy_603", (0.9985390234525182, 0.9547097270280661, 0.6803537101114956))
cmd.color("dummy_603", "test2 and state 604", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 604")

cmd.set_color("dummy_604", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_604", "test2 and state 605", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 605")

cmd.set_color("dummy_605", (0.9966935793925413, 0.8975009611687812, 0.5936178392925797))
cmd.color("dummy_605", "test2 and state 606", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 606")

cmd.set_color("dummy_606", (0.9301038062283737, 0.3711649365628604, 0.23690888119953865))
cmd.color("dummy_606", "test2 and state 607", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 607")

cmd.set_color("dummy_607", (0.9390234525182622, 0.3899269511726259, 0.245520953479431))
cmd.color("dummy_607", "test2 and state 608", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 608")

cmd.set_color("dummy_608", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_608", "test2 and state 609", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 609")

cmd.set_color("dummy_609", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_609", "test2 and state 610", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 610")

cmd.set_color("dummy_610", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_610", "test2 and state 611", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 611")

cmd.set_color("dummy_611", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_611", "test2 and state 612", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 612")

cmd.set_color("dummy_612", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_612", "test2 and state 613", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 613")

cmd.set_color("dummy_613", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_613", "test2 and state 614", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 614")

cmd.set_color("dummy_614", (0.9933102652825836, 0.7400230680507497, 0.43460207612456747))
cmd.color("dummy_614", "test2 and state 615", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 615")

cmd.set_color("dummy_615", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_615", "test2 and state 616", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 616")

cmd.set_color("dummy_616", (0.9979238754325259, 0.9356401384083045, 0.6514417531718569))
cmd.color("dummy_616", "test2 and state 617", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 617")

cmd.set_color("dummy_617", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_617", "test2 and state 618", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 618")

cmd.set_color("dummy_618", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_618", "test2 and state 619", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 619")

cmd.set_color("dummy_619", (0.996078431372549, 0.8784313725490196, 0.5647058823529412))
cmd.color("dummy_619", "test2 and state 620", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 620")

cmd.set_color("dummy_620", (0.9988465974625145, 0.9642445213379469, 0.6948096885813149))
cmd.color("dummy_620", "test2 and state 621", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 621")

cmd.set_color("dummy_621", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_621", "test2 and state 622", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 622")

cmd.set_color("dummy_622", (0.994079200307574, 0.7784698193002692, 0.4707420222991157))
cmd.color("dummy_622", "test2 and state 623", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 623")

cmd.set_color("dummy_623", (0.9953094963475586, 0.8399846212995001, 0.5285659361783929))
cmd.color("dummy_623", "test2 and state 624", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 624")

cmd.set_color("dummy_624", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_624", "test2 and state 625", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 625")

cmd.set_color("dummy_625", (0.9982314494425221, 0.9451749327181853, 0.6658977316416763))
cmd.color("dummy_625", "test2 and state 626", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 626")

cmd.set_color("dummy_626", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_626", "test2 and state 627", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 627")

cmd.set_color("dummy_627", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_627", "test2 and state 628", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 628")

cmd.set_color("dummy_628", (0.3031141868512111, 0.49757785467128035, 0.7259515570934256))
cmd.color("dummy_628", "test2 and state 629", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 629")

cmd.set_color("dummy_629", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_629", "test2 and state 630", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 630")

cmd.set_color("dummy_630", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_630", "test2 and state 631", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 631")

cmd.set_color("dummy_631", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_631", "test2 and state 632", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 632")

cmd.set_color("dummy_632", (0.19215686274509805, 0.21176470588235294, 0.5843137254901961))
cmd.color("dummy_632", "test2 and state 633", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 633")

cmd.set_color("dummy_633", (0.9950019223375625, 0.8246059207996924, 0.5141099577085736))
cmd.color("dummy_633", "test2 and state 634", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 634")

cmd.set_color("dummy_634", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_634", "test2 and state 635", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 635")

cmd.set_color("dummy_635", (0.4187620146097656, 0.635371011149558, 0.7973087274125337))
cmd.color("dummy_635", "test2 and state 636", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 636")

cmd.set_color("dummy_636", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_636", "test2 and state 637", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 637")

cmd.set_color("dummy_637", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_637", "test2 and state 638", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 638")

cmd.set_color("dummy_638", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_638", "test2 and state 639", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 639")

cmd.set_color("dummy_639", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_639", "test2 and state 640", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 640")

cmd.set_color("dummy_640", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_640", "test2 and state 641", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 641")

cmd.set_color("dummy_641", (0.9997693194925029, 0.9928489042675894, 0.7381776239907728))
cmd.color("dummy_641", "test2 and state 642", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 642")

cmd.set_color("dummy_642", (0.9945405613225683, 0.8015378700499808, 0.4924259900038447))
cmd.color("dummy_642", "test2 and state 643", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 643")

cmd.set_color("dummy_643", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_643", "test2 and state 644", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 644")

cmd.set_color("dummy_644", (0.9936178392925799, 0.7554017685505575, 0.44905805459438675))
cmd.color("dummy_644", "test2 and state 645", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 645")

cmd.set_color("dummy_645", (0.5310265282583623, 0.7393310265282584, 0.8528258362168397))
cmd.color("dummy_645", "test2 and state 646", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 646")

cmd.set_color("dummy_646", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_646", "test2 and state 647", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 647")

cmd.set_color("dummy_647", (0.9922337562475971, 0.6861976163014225, 0.3840061514801999))
cmd.color("dummy_647", "test2 and state 648", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 648")

cmd.set_color("dummy_648", (0.9547097270280662, 0.9824682814302191, 0.8322952710495962))
cmd.color("dummy_648", "test2 and state 649", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 649")

cmd.set_color("dummy_649", (0.9982314494425221, 0.9451749327181853, 0.6658977316416763))
cmd.color("dummy_649", "test2 and state 650", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 650")

cmd.set_color("dummy_650", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_650", "test2 and state 651", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 651")

cmd.set_color("dummy_651", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_651", "test2 and state 652", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 652")

cmd.set_color("dummy_652", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_652", "test2 and state 653", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 653")

cmd.set_color("dummy_653", (0.9966935793925413, 0.8975009611687812, 0.5936178392925797))
cmd.color("dummy_653", "test2 and state 654", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 654")

cmd.set_color("dummy_654", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_654", "test2 and state 655", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 655")

cmd.set_color("dummy_655", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_655", "test2 and state 656", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 656")

cmd.set_color("dummy_656", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_656", "test2 and state 657", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 657")

cmd.set_color("dummy_657", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_657", "test2 and state 658", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 658")

cmd.set_color("dummy_658", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_658", "test2 and state 659", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 659")

cmd.set_color("dummy_659", (0.996078431372549, 0.8784313725490196, 0.5647058823529412))
cmd.color("dummy_659", "test2 and state 660", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 660")

cmd.set_color("dummy_660", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_660", "test2 and state 661", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 661")

cmd.set_color("dummy_661", (0.5225682429834679, 0.7325643983083431, 0.8491349480968858))
cmd.color("dummy_661", "test2 and state 662", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 662")

cmd.set_color("dummy_662", (0.2742022299115725, 0.4631295655517109, 0.7081122645136487))
cmd.color("dummy_662", "test2 and state 663", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 663")

cmd.set_color("dummy_663", (0.9982314494425221, 0.9451749327181853, 0.6658977316416763))
cmd.color("dummy_663", "test2 and state 664", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 664")

cmd.set_color("dummy_664", (0.9988465974625145, 0.9642445213379469, 0.6948096885813149))
cmd.color("dummy_664", "test2 and state 665", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 665")

cmd.set_color("dummy_665", (0.6828143021914649, 0.8569780853517878, 0.9171856978085352))
cmd.color("dummy_665", "test2 and state 666", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 666")

cmd.set_color("dummy_666", (0.5141099577085737, 0.7257977700884276, 0.845444059976932))
cmd.color("dummy_666", "test2 and state 667", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 667")

cmd.set_color("dummy_667", (0.9990003844675125, 0.9690119184928874, 0.7020376778162245))
cmd.color("dummy_667", "test2 and state 668", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 668")

cmd.set_color("dummy_668", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_668", "test2 and state 669", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 669")

cmd.set_color("dummy_669", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_669", "test2 and state 670", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 670")

cmd.set_color("dummy_670", (0.9936178392925799, 0.7554017685505575, 0.44905805459438675))
cmd.color("dummy_670", "test2 and state 671", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 671")

cmd.set_color("dummy_671", (0.9950019223375625, 0.8246059207996924, 0.5141099577085736))
cmd.color("dummy_671", "test2 and state 672", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 672")

cmd.set_color("dummy_672", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_672", "test2 and state 673", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 673")

cmd.set_color("dummy_673", (0.9707035755478662, 0.5274125336409073, 0.308881199538639))
cmd.color("dummy_673", "test2 and state 674", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 674")

cmd.set_color("dummy_674", (0.9790080738177624, 0.5873894655901576, 0.336562860438293))
cmd.color("dummy_674", "test2 and state 675", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 675")

cmd.set_color("dummy_675", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_675", "test2 and state 676", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 676")

cmd.set_color("dummy_676", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_676", "test2 and state 677", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 677")

cmd.set_color("dummy_677", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_677", "test2 and state 678", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 678")

cmd.set_color("dummy_678", (0.9951557093425606, 0.8322952710495963, 0.5213379469434832))
cmd.color("dummy_678", "test2 and state 679", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 679")

cmd.set_color("dummy_679", (0.9434832756632064, 0.39930795847750866, 0.24982698961937716))
cmd.color("dummy_679", "test2 and state 680", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 680")

cmd.set_color("dummy_680", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_680", "test2 and state 681", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 681")

cmd.set_color("dummy_681", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_681", "test2 and state 682", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 682")

cmd.set_color("dummy_682", (0.5986928104575164, 0.7934640522875818, 0.8823529411764706))
cmd.color("dummy_682", "test2 and state 683", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 683")

cmd.set_color("dummy_683", (0.998077662437524, 0.9404075355632449, 0.6586697424067667))
cmd.color("dummy_683", "test2 and state 684", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 684")

cmd.set_color("dummy_684", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_684", "test2 and state 685", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 685")

cmd.set_color("dummy_685", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_685", "test2 and state 686", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 686")

cmd.set_color("dummy_686", (0.9547097270280662, 0.9824682814302191, 0.8322952710495962))
cmd.color("dummy_686", "test2 and state 687", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 687")

cmd.set_color("dummy_687", (0.9886966551326413, 0.657362552864283, 0.36885813148788926))
cmd.color("dummy_687", "test2 and state 688", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 688")

cmd.set_color("dummy_688", (0.9451749327181853, 0.9787773933102653, 0.8498269896193771))
cmd.color("dummy_688", "test2 and state 689", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 689")

cmd.set_color("dummy_689", (0.8587466359092657, 0.22106881968473663, 0.16801230296039987))
cmd.color("dummy_689", "test2 and state 690", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 690")

cmd.set_color("dummy_690", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_690", "test2 and state 691", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 691")

cmd.set_color("dummy_691", (0.4971933871587852, 0.7122645136485968, 0.8380622837370243))
cmd.color("dummy_691", "test2 and state 692", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 692")

cmd.set_color("dummy_692", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_692", "test2 and state 693", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 693")

cmd.set_color("dummy_693", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_693", "test2 and state 694", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 694")

cmd.set_color("dummy_694", (0.7724721261053442, 0.900961168781238, 0.942560553633218))
cmd.color("dummy_694", "test2 and state 695", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 695")

cmd.set_color("dummy_695", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_695", "test2 and state 696", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 696")

cmd.set_color("dummy_696", (0.9922337562475971, 0.6861976163014225, 0.3840061514801999))
cmd.color("dummy_696", "test2 and state 697", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 697")

cmd.set_color("dummy_697", (0.8539792387543255, 0.9409457900807382, 0.9656286043829296))
cmd.color("dummy_697", "test2 and state 698", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 698")

cmd.set_color("dummy_698", (0.4259900038446752, 0.6439830834294503, 0.801768550557478))
cmd.color("dummy_698", "test2 and state 699", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 699")

cmd.set_color("dummy_699", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_699", "test2 and state 700", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 700")

cmd.set_color("dummy_700", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_700", "test2 and state 701", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 701")

cmd.set_color("dummy_701", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_701", "test2 and state 702", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 702")

cmd.set_color("dummy_702", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_702", "test2 and state 703", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 703")

cmd.set_color("dummy_703", (0.8539792387543255, 0.9409457900807382, 0.9656286043829296))
cmd.color("dummy_703", "test2 and state 704", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 704")

cmd.set_color("dummy_704", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_704", "test2 and state 705", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 705")

cmd.set_color("dummy_705", (0.9833141099577086, 0.9935409457900808, 0.7797001153402537))
cmd.color("dummy_705", "test2 and state 706", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 706")

cmd.set_color("dummy_706", (0.31757016532103044, 0.5148019992310651, 0.7348712033833141))
cmd.color("dummy_706", "test2 and state 707", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 707")

cmd.set_color("dummy_707", (0.9983852364475202, 0.9499423298731258, 0.673125720876586))
cmd.color("dummy_707", "test2 and state 708", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 708")

cmd.set_color("dummy_708", (0.9973087274125336, 0.9165705497885428, 0.6225297962322184))
cmd.color("dummy_708", "test2 and state 709", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 709")

cmd.set_color("dummy_709", (0.9928489042675894, 0.9972318339100346, 0.7621683967704729))
cmd.color("dummy_709", "test2 and state 710", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 710")

cmd.set_color("dummy_710", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_710", "test2 and state 711", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 711")

cmd.set_color("dummy_711", (0.2413687043444829, 0.3667820069204153, 0.6605920799692426))
cmd.color("dummy_711", "test2 and state 712", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 712")

cmd.set_color("dummy_712", (0.9665513264129182, 0.4974240676662822, 0.295040369088812))
cmd.color("dummy_712", "test2 and state 713", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 713")

cmd.set_color("dummy_713", (0.7162629757785467, 0.06643598615916954, 0.15040369088811995))
cmd.color("dummy_713", "test2 and state 714", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 714")

cmd.set_color("dummy_714", (0.5225682429834679, 0.7325643983083431, 0.8491349480968858))
cmd.color("dummy_714", "test2 and state 715", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 715")

cmd.set_color("dummy_715", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_715", "test2 and state 716", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 716")

cmd.set_color("dummy_716", (0.9986928104575163, 0.9594771241830066, 0.6875816993464052))
cmd.color("dummy_716", "test2 and state 717", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 717")

cmd.set_color("dummy_717", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_717", "test2 and state 718", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 718")

cmd.set_color("dummy_718", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_718", "test2 and state 719", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 719")

cmd.set_color("dummy_719", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_719", "test2 and state 720", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 720")

cmd.set_color("dummy_720", (0.9679354094579009, 0.5074202229911575, 0.2996539792387545))
cmd.color("dummy_720", "test2 and state 721", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 721")

cmd.set_color("dummy_721", (0.9923875432525952, 0.6938869665513264, 0.39123414071510954))
cmd.color("dummy_721", "test2 and state 722", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 722")

cmd.set_color("dummy_722", (0.9748558246828143, 0.5574009996155325, 0.322722029988466))
cmd.color("dummy_722", "test2 and state 723", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 723")

cmd.set_color("dummy_723", (0.9974625144175318, 0.9213379469434833, 0.629757785467128))
cmd.color("dummy_723", "test2 and state 724", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 724")

cmd.set_color("dummy_724", (0.9933102652825836, 0.7400230680507497, 0.43460207612456747))
cmd.color("dummy_724", "test2 and state 725", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 725")

cmd.set_color("dummy_725", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_725", "test2 and state 726", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 726")

cmd.set_color("dummy_726", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_726", "test2 and state 727", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 727")

cmd.set_color("dummy_727", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_727", "test2 and state 728", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 728")

cmd.set_color("dummy_728", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_728", "test2 and state 729", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 729")

cmd.set_color("dummy_729", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_729", "test2 and state 730", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 730")

cmd.set_color("dummy_730", (0.4971933871587852, 0.7122645136485968, 0.8380622837370243))
cmd.color("dummy_730", "test2 and state 731", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 731")

cmd.set_color("dummy_731", (0.9845444059976932, 0.6273740868896578, 0.3550173010380623))
cmd.color("dummy_731", "test2 and state 732", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 732")

cmd.set_color("dummy_732", (0.9962322183775472, 0.88319876970396, 0.5719338715878508))
cmd.color("dummy_732", "test2 and state 733", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 733")

cmd.set_color("dummy_733", (0.9928489042675894, 0.9972318339100346, 0.7621683967704729))
cmd.color("dummy_733", "test2 and state 734", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 734")

cmd.set_color("dummy_734", (0.9993079584775086, 0.9785467128027683, 0.716493656286044))
cmd.color("dummy_734", "test2 and state 735", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 735")

cmd.set_color("dummy_735", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_735", "test2 and state 736", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 736")

cmd.set_color("dummy_736", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_736", "test2 and state 737", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 737")

cmd.set_color("dummy_737", (0.9213379469434834, 0.9695501730103806, 0.8936562860438291))
cmd.color("dummy_737", "test2 and state 738", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 738")

cmd.set_color("dummy_738", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_738", "test2 and state 739", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 739")

cmd.set_color("dummy_739", (0.43321799307958486, 0.6525951557093427, 0.8062283737024222))
cmd.color("dummy_739", "test2 and state 740", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 740")

cmd.set_color("dummy_740", (0.998077662437524, 0.9404075355632449, 0.6586697424067667))
cmd.color("dummy_740", "test2 and state 741", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 741")

cmd.set_color("dummy_741", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_741", "test2 and state 742", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 742")

cmd.set_color("dummy_742", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_742", "test2 and state 743", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 743")

cmd.set_color("dummy_743", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_743", "test2 and state 744", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 744")

cmd.set_color("dummy_744", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_744", "test2 and state 745", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 745")

cmd.set_color("dummy_745", (0.8008458285274894, 0.14763552479815456, 0.1520953479430988))
cmd.color("dummy_745", "test2 and state 746", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 746")

cmd.set_color("dummy_746", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_746", "test2 and state 747", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 747")

cmd.set_color("dummy_747", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_747", "test2 and state 748", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 748")

cmd.set_color("dummy_748", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_748", "test2 and state 749", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 749")

cmd.set_color("dummy_749", (0.4633602460592081, 0.6851980007689351, 0.8232987312572088))
cmd.color("dummy_749", "test2 and state 750", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 750")

cmd.set_color("dummy_750", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_750", "test2 and state 751", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 751")

cmd.set_color("dummy_751", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_751", "test2 and state 752", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 752")

cmd.set_color("dummy_752", (0.9976163014225298, 0.9261053440984237, 0.6369857747020377))
cmd.color("dummy_752", "test2 and state 753", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 753")

cmd.set_color("dummy_753", (0.9499423298731258, 0.9806228373702423, 0.8410611303344866))
cmd.color("dummy_753", "test2 and state 754", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 754")

cmd.set_color("dummy_754", (0.9900807381776241, 0.6673587081891583, 0.3734717416378317))
cmd.color("dummy_754", "test2 and state 755", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 755")

cmd.set_color("dummy_755", (0.996078431372549, 0.8784313725490196, 0.5647058823529412))
cmd.color("dummy_755", "test2 and state 756", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 756")

cmd.set_color("dummy_756", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_756", "test2 and state 757", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 757")

cmd.set_color("dummy_757", (0.5564013840830451, 0.7596309111880047, 0.8638985005767013))
cmd.color("dummy_757", "test2 and state 758", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 758")

cmd.set_color("dummy_758", (0.9988465974625145, 0.9642445213379469, 0.6948096885813149))
cmd.color("dummy_758", "test2 and state 759", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 759")

cmd.set_color("dummy_759", (0.9974625144175318, 0.9213379469434833, 0.629757785467128))
cmd.color("dummy_759", "test2 and state 760", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 760")

cmd.set_color("dummy_760", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_760", "test2 and state 761", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 761")

cmd.set_color("dummy_761", (0.3898500576701269, 0.6009227220299886, 0.7794694348327567))
cmd.color("dummy_761", "test2 and state 762", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 762")

cmd.set_color("dummy_762", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_762", "test2 and state 763", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 763")

cmd.set_color("dummy_763", (0.9845444059976932, 0.6273740868896578, 0.3550173010380623))
cmd.color("dummy_763", "test2 and state 764", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 764")

cmd.set_color("dummy_764", (0.9966935793925413, 0.8975009611687812, 0.5936178392925797))
cmd.color("dummy_764", "test2 and state 765", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 765")

cmd.set_color("dummy_765", (0.8702806612841217, 0.9489427143406383, 0.9702422145328721))
cmd.color("dummy_765", "test2 and state 766", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 766")

cmd.set_color("dummy_766", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_766", "test2 and state 767", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 767")

cmd.set_color("dummy_767", (0.4187620146097656, 0.635371011149558, 0.7973087274125337))
cmd.color("dummy_767", "test2 and state 768", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 768")

cmd.set_color("dummy_768", (0.38262206843521723, 0.5923106497500962, 0.7750096116878125))
cmd.color("dummy_768", "test2 and state 769", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 769")

cmd.set_color("dummy_769", (0.9951557093425606, 0.8322952710495963, 0.5213379469434832))
cmd.color("dummy_769", "test2 and state 770", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 770")

cmd.set_color("dummy_770", (0.8702806612841217, 0.9489427143406383, 0.9702422145328721))
cmd.color("dummy_770", "test2 and state 771", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 771")

cmd.set_color("dummy_771", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_771", "test2 and state 772", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 772")

cmd.set_color("dummy_772", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_772", "test2 and state 773", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 773")

cmd.set_color("dummy_773", (0.9859284890426759, 0.6373702422145329, 0.3596309111880046))
cmd.color("dummy_773", "test2 and state 774", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 774")

cmd.set_color("dummy_774", (0.8453671664744329, 0.19292579777008845, 0.15509419454056134))
cmd.color("dummy_774", "test2 and state 775", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 775")

cmd.set_color("dummy_775", (0.9122645136485967, 0.33364090734332946, 0.21968473663975396))
cmd.color("dummy_775", "test2 and state 776", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 776")

cmd.set_color("dummy_776", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_776", "test2 and state 777", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 777")

cmd.set_color("dummy_777", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_777", "test2 and state 778", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 778")

cmd.set_color("dummy_778", (0.47181853133410234, 0.6919646289888505, 0.8269896193771626))
cmd.color("dummy_778", "test2 and state 779", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 779")

cmd.set_color("dummy_779", (0.9873125720876587, 0.647366397539408, 0.36424452133794694))
cmd.color("dummy_779", "test2 and state 780", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 780")

cmd.set_color("dummy_780", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_780", "test2 and state 781", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 781")

cmd.set_color("dummy_781", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_781", "test2 and state 782", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 782")

cmd.set_color("dummy_782", (0.5986928104575164, 0.7934640522875818, 0.8823529411764706))
cmd.color("dummy_782", "test2 and state 783", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 783")

cmd.set_color("dummy_783", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_783", "test2 and state 784", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 784")

cmd.set_color("dummy_784", (0.8498269896193773, 0.20230680507497142, 0.15940023068050763))
cmd.color("dummy_784", "test2 and state 785", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 785")

cmd.set_color("dummy_785", (0.9211841599384853, 0.3524029219530952, 0.22829680891964643))
cmd.color("dummy_785", "test2 and state 786", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 786")

cmd.set_color("dummy_786", (0.9923875432525952, 0.6938869665513264, 0.39123414071510954))
cmd.color("dummy_786", "test2 and state 787", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 787")

cmd.set_color("dummy_787", (0.39707804690503656, 0.6095347943098809, 0.783929257977701))
cmd.color("dummy_787", "test2 and state 788", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 788")

cmd.set_color("dummy_788", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_788", "test2 and state 789", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 789")

cmd.set_color("dummy_789", (0.5310265282583623, 0.7393310265282584, 0.8528258362168397))
cmd.color("dummy_789", "test2 and state 790", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 790")

cmd.set_color("dummy_790", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_790", "test2 and state 791", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 791")

cmd.set_color("dummy_791", (0.994079200307574, 0.7784698193002692, 0.4707420222991157))
cmd.color("dummy_791", "test2 and state 792", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 792")

cmd.set_color("dummy_792", (0.998077662437524, 0.9404075355632449, 0.6586697424067667))
cmd.color("dummy_792", "test2 and state 793", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 793")

cmd.set_color("dummy_793", (0.9997693194925029, 0.9928489042675894, 0.7381776239907728))
cmd.color("dummy_793", "test2 and state 794", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 794")

cmd.set_color("dummy_794", (0.9973087274125336, 0.9165705497885428, 0.6225297962322184))
cmd.color("dummy_794", "test2 and state 795", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 795")

cmd.set_color("dummy_795", (0.9933102652825836, 0.7400230680507497, 0.43460207612456747))
cmd.color("dummy_795", "test2 and state 796", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 796")

cmd.set_color("dummy_796", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_796", "test2 and state 797", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 797")

cmd.set_color("dummy_797", (0.5141099577085737, 0.7257977700884276, 0.845444059976932))
cmd.color("dummy_797", "test2 and state 798", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 798")

cmd.set_color("dummy_798", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_798", "test2 and state 799", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 799")

cmd.set_color("dummy_799", (0.9434832756632064, 0.39930795847750866, 0.24982698961937716))
cmd.color("dummy_799", "test2 and state 800", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 800")

cmd.set_color("dummy_800", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_800", "test2 and state 801", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 801")

cmd.set_color("dummy_801", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_801", "test2 and state 802", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 802")

cmd.set_color("dummy_802", (0.8988850442137639, 0.3054978854286813, 0.20676662821991543))
cmd.color("dummy_802", "test2 and state 803", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 803")

cmd.set_color("dummy_803", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_803", "test2 and state 804", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 804")

cmd.set_color("dummy_804", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_804", "test2 and state 805", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 805")

cmd.set_color("dummy_805", (0.7072664359861592, 0.8689734717416379, 0.9241061130334487))
cmd.color("dummy_805", "test2 and state 806", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 806")

cmd.set_color("dummy_806", (0.5479430988081508, 0.7528642829680893, 0.8602076124567475))
cmd.color("dummy_806", "test2 and state 807", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 807")

cmd.set_color("dummy_807", (0.9976163014225298, 0.9261053440984237, 0.6369857747020377))
cmd.color("dummy_807", "test2 and state 808", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 808")

cmd.set_color("dummy_808", (0.9833141099577086, 0.9935409457900808, 0.7797001153402537))
cmd.color("dummy_808", "test2 and state 809", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 809")

cmd.set_color("dummy_809", (0.9979238754325259, 0.9356401384083045, 0.6514417531718569))
cmd.color("dummy_809", "test2 and state 810", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 810")

cmd.set_color("dummy_810", (0.9925413302575933, 0.7015763168012303, 0.3984621299500192))
cmd.color("dummy_810", "test2 and state 811", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 811")

cmd.set_color("dummy_811", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_811", "test2 and state 812", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 812")

cmd.set_color("dummy_812", (0.996078431372549, 0.8784313725490196, 0.5647058823529412))
cmd.color("dummy_812", "test2 and state 813", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 813")

cmd.set_color("dummy_813", (0.4802768166089966, 0.6987312572087659, 0.8306805074971165))
cmd.color("dummy_813", "test2 and state 814", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 814")

cmd.set_color("dummy_814", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_814", "test2 and state 815", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 815")

cmd.set_color("dummy_815", (0.5056516724336794, 0.7190311418685121, 0.8417531718569781))
cmd.color("dummy_815", "test2 and state 816", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 816")

cmd.set_color("dummy_816", (0.9308727412533642, 0.9732410611303345, 0.8761245674740483))
cmd.color("dummy_816", "test2 and state 817", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 817")

cmd.set_color("dummy_817", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_817", "test2 and state 818", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 818")

cmd.set_color("dummy_818", (0.31034217608612075, 0.5061899269511727, 0.7304113802383699))
cmd.color("dummy_818", "test2 and state 819", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 819")

cmd.set_color("dummy_819", (0.43321799307958486, 0.6525951557093427, 0.8062283737024222))
cmd.color("dummy_819", "test2 and state 820", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 820")

cmd.set_color("dummy_820", (0.9999231064975009, 0.9976163014225298, 0.7454056132256824))
cmd.color("dummy_820", "test2 and state 821", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 821")

cmd.set_color("dummy_821", (0.43321799307958486, 0.6525951557093427, 0.8062283737024222))
cmd.color("dummy_821", "test2 and state 822", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 822")

cmd.set_color("dummy_822", (0.9547097270280662, 0.9824682814302191, 0.8322952710495962))
cmd.color("dummy_822", "test2 and state 823", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 823")

cmd.set_color("dummy_823", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_823", "test2 and state 824", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 824")

cmd.set_color("dummy_824", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_824", "test2 and state 825", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 825")

cmd.set_color("dummy_825", (0.3898500576701269, 0.6009227220299886, 0.7794694348327567))
cmd.color("dummy_825", "test2 and state 826", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 826")

cmd.set_color("dummy_826", (0.8831987697039602, 0.9547866205305652, 0.9637831603229524))
cmd.color("dummy_826", "test2 and state 827", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 827")

cmd.set_color("dummy_827", (0.9979238754325259, 0.9356401384083045, 0.6514417531718569))
cmd.color("dummy_827", "test2 and state 828", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 828")

cmd.set_color("dummy_828", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_828", "test2 and state 829", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 829")

cmd.set_color("dummy_829", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_829", "test2 and state 830", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 830")

cmd.set_color("dummy_830", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_830", "test2 and state 831", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 831")

cmd.set_color("dummy_831", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_831", "test2 and state 832", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 832")

cmd.set_color("dummy_832", (0.4887351018838909, 0.7054978854286814, 0.8343713956170704))
cmd.color("dummy_832", "test2 and state 833", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 833")

cmd.set_color("dummy_833", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_833", "test2 and state 834", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 834")

cmd.set_color("dummy_834", (0.4802768166089966, 0.6987312572087659, 0.8306805074971165))
cmd.color("dummy_834", "test2 and state 835", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 835")

cmd.set_color("dummy_835", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_835", "test2 and state 836", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 836")

cmd.set_color("dummy_836", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_836", "test2 and state 837", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 837")

cmd.set_color("dummy_837", (0.9950019223375625, 0.8246059207996924, 0.5141099577085736))
cmd.color("dummy_837", "test2 and state 838", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 838")

cmd.set_color("dummy_838", (0.952402921953095, 0.4180699730872741, 0.2584390618992695))
cmd.color("dummy_838", "test2 and state 839", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 839")

cmd.set_color("dummy_839", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_839", "test2 and state 840", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 840")

cmd.set_color("dummy_840", (0.6240676662821992, 0.8137639369473281, 0.8934256055363322))
cmd.color("dummy_840", "test2 and state 841", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 841")

cmd.set_color("dummy_841", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_841", "test2 and state 842", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 842")

cmd.set_color("dummy_842", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_842", "test2 and state 843", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 843")

cmd.set_color("dummy_843", (0.6240676662821992, 0.8137639369473281, 0.8934256055363322))
cmd.color("dummy_843", "test2 and state 844", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 844")

cmd.set_color("dummy_844", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_844", "test2 and state 845", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 845")

cmd.set_color("dummy_845", (0.9928489042675894, 0.7169550173010381, 0.4129181084198385))
cmd.color("dummy_845", "test2 and state 846", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 846")

cmd.set_color("dummy_846", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_846", "test2 and state 847", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 847")

cmd.set_color("dummy_847", (0.9954632833525567, 0.847673971549404, 0.5357939254133026))
cmd.color("dummy_847", "test2 and state 848", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 848")

cmd.set_color("dummy_848", (0.96239907727797, 0.46743560169165704, 0.281199538638985))
cmd.color("dummy_848", "test2 and state 849", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 849")

cmd.set_color("dummy_849", (0.9167243367935409, 0.3430219146482122, 0.22399077277970011))
cmd.color("dummy_849", "test2 and state 850", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 850")

cmd.set_color("dummy_850", (0.4404459823144945, 0.661207227989235, 0.8106881968473664))
cmd.color("dummy_850", "test2 and state 851", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 851")

cmd.set_color("dummy_851", (0.8944252210688197, 0.2961168781237985, 0.20246059207996925))
cmd.color("dummy_851", "test2 and state 852", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 852")

cmd.set_color("dummy_852", (0.756170703575548, 0.8929642445213379, 0.9379469434832757))
cmd.color("dummy_852", "test2 and state 853", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 853")

cmd.set_color("dummy_853", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_853", "test2 and state 854", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 854")

cmd.set_color("dummy_854", (0.41153402537485584, 0.6267589388696656, 0.7928489042675894))
cmd.color("dummy_854", "test2 and state 855", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 855")

cmd.set_color("dummy_855", (0.9568627450980393, 0.42745098039215684, 0.2627450980392157))
cmd.color("dummy_855", "test2 and state 856", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 856")

cmd.set_color("dummy_856", (0.9930026912725874, 0.724644367550942, 0.4201460976547482))
cmd.color("dummy_856", "test2 and state 857", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 857")

cmd.set_color("dummy_857", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_857", "test2 and state 858", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 858")

cmd.set_color("dummy_858", (0.45490196078431383, 0.6784313725490198, 0.8196078431372549))
cmd.color("dummy_858", "test2 and state 859", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 859")

cmd.set_color("dummy_859", (0.9982314494425221, 0.9451749327181853, 0.6658977316416763))
cmd.color("dummy_859", "test2 and state 860", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 860")

cmd.set_color("dummy_860", (0.9665513264129182, 0.4974240676662822, 0.295040369088812))
cmd.color("dummy_860", "test2 and state 861", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 861")

cmd.set_color("dummy_861", (0.4259900038446752, 0.6439830834294503, 0.801768550557478))
cmd.color("dummy_861", "test2 and state 862", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 862")

cmd.set_color("dummy_862", (0.9451749327181853, 0.9787773933102653, 0.8498269896193771))
cmd.color("dummy_862", "test2 and state 863", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 863")

cmd.set_color("dummy_863", (0.9988465974625145, 0.9642445213379469, 0.6948096885813149))
cmd.color("dummy_863", "test2 and state 864", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 864")

cmd.set_color("dummy_864", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_864", "test2 and state 865", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 865")

cmd.set_color("dummy_865", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_865", "test2 and state 866", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 866")

cmd.set_color("dummy_866", (0.8702806612841217, 0.9489427143406383, 0.9702422145328721))
cmd.color("dummy_866", "test2 and state 867", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 867")

cmd.set_color("dummy_867", (0.9983852364475202, 0.9499423298731258, 0.673125720876586))
cmd.color("dummy_867", "test2 and state 868", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 868")

cmd.set_color("dummy_868", (0.4259900038446752, 0.6439830834294503, 0.801768550557478))
cmd.color("dummy_868", "test2 and state 869", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 869")

cmd.set_color("dummy_869", (0.21983852364475204, 0.298961937716263, 0.6272202998846598))
cmd.color("dummy_869", "test2 and state 870", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 870")

cmd.set_color("dummy_870", (0.9976163014225298, 0.9990772779700116, 0.7534025374855824))
cmd.color("dummy_870", "test2 and state 871", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 871")

cmd.set_color("dummy_871", (0.9451749327181853, 0.9787773933102653, 0.8498269896193771))
cmd.color("dummy_871", "test2 and state 872", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 872")

cmd.set_color("dummy_872", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_872", "test2 and state 873", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 873")

cmd.set_color("dummy_873", (0.5141099577085737, 0.7257977700884276, 0.845444059976932))
cmd.color("dummy_873", "test2 and state 874", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 874")

cmd.set_color("dummy_874", (0.38262206843521723, 0.5923106497500962, 0.7750096116878125))
cmd.color("dummy_874", "test2 and state 875", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 875")

cmd.set_color("dummy_875", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_875", "test2 and state 876", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 876")

cmd.set_color("dummy_876", (0.34648212226066905, 0.5492502883506345, 0.7527104959630911))
cmd.color("dummy_876", "test2 and state 877", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 877")

cmd.set_color("dummy_877", (0.5564013840830451, 0.7596309111880047, 0.8638985005767013))
cmd.color("dummy_877", "test2 and state 878", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 878")

cmd.set_color("dummy_878", (0.9936178392925799, 0.7554017685505575, 0.44905805459438675))
cmd.color("dummy_878", "test2 and state 879", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 879")

cmd.set_color("dummy_879", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_879", "test2 and state 880", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 880")

cmd.set_color("dummy_880", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_880", "test2 and state 881", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 881")

cmd.set_color("dummy_881", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_881", "test2 and state 882", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 882")

cmd.set_color("dummy_882", (0.9928489042675894, 0.9972318339100346, 0.7621683967704729))
cmd.color("dummy_882", "test2 and state 883", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 883")

cmd.set_color("dummy_883", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_883", "test2 and state 884", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 884")

cmd.set_color("dummy_884", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_884", "test2 and state 885", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 885")

cmd.set_color("dummy_885", (0.9833141099577086, 0.9935409457900808, 0.7797001153402537))
cmd.color("dummy_885", "test2 and state 886", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 886")

cmd.set_color("dummy_886", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_886", "test2 and state 887", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 887")

cmd.set_color("dummy_887", (0.9985390234525182, 0.9547097270280661, 0.6803537101114956))
cmd.color("dummy_887", "test2 and state 888", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 888")

cmd.set_color("dummy_888", (0.9994617454825068, 0.9833141099577085, 0.7237216455209535))
cmd.color("dummy_888", "test2 and state 889", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 889")

cmd.set_color("dummy_889", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_889", "test2 and state 890", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 890")

cmd.set_color("dummy_890", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_890", "test2 and state 891", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 891")

cmd.set_color("dummy_891", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_891", "test2 and state 892", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 892")

cmd.set_color("dummy_892", (0.940407535563245, 0.9769319492502884, 0.8585928489042675))
cmd.color("dummy_892", "test2 and state 893", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 893")

cmd.set_color("dummy_893", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_893", "test2 and state 894", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 894")

cmd.set_color("dummy_894", (0.8702806612841217, 0.9489427143406383, 0.9702422145328721))
cmd.color("dummy_894", "test2 and state 895", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 895")

cmd.set_color("dummy_895", (0.9070357554786621, 0.9640138408304498, 0.9199538638985004))
cmd.color("dummy_895", "test2 and state 896", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 896")

cmd.set_color("dummy_896", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_896", "test2 and state 897", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 897")

cmd.set_color("dummy_897", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_897", "test2 and state 898", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 898")

cmd.set_color("dummy_898", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_898", "test2 and state 899", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 899")

cmd.set_color("dummy_899", (0.39707804690503656, 0.6095347943098809, 0.783929257977701))
cmd.color("dummy_899", "test2 and state 900", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 900")

cmd.set_color("dummy_900", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_900", "test2 and state 901", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 901")

cmd.set_color("dummy_901", (0.31034217608612075, 0.5061899269511727, 0.7304113802383699))
cmd.color("dummy_901", "test2 and state 902", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 902")

cmd.set_color("dummy_902", (0.5141099577085737, 0.7257977700884276, 0.845444059976932))
cmd.color("dummy_902", "test2 and state 903", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 903")

cmd.set_color("dummy_903", (0.7643214148404461, 0.896962706651288, 0.9402537485582468))
cmd.color("dummy_903", "test2 and state 904", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 904")

cmd.set_color("dummy_904", (0.9803921568627452, 0.5973856209150327, 0.3411764705882353))
cmd.color("dummy_904", "test2 and state 905", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 905")

cmd.set_color("dummy_905", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_905", "test2 and state 906", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 906")

cmd.set_color("dummy_906", (0.9451749327181853, 0.9787773933102653, 0.8498269896193771))
cmd.color("dummy_906", "test2 and state 907", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 907")

cmd.set_color("dummy_907", (0.9922337562475971, 0.6861976163014225, 0.3840061514801999))
cmd.color("dummy_907", "test2 and state 908", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 908")

cmd.set_color("dummy_908", (0.8008458285274894, 0.14763552479815456, 0.1520953479430988))
cmd.color("dummy_908", "test2 and state 909", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 909")

cmd.set_color("dummy_909", (0.952402921953095, 0.4180699730872741, 0.2584390618992695))
cmd.color("dummy_909", "test2 and state 910", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 910")

cmd.set_color("dummy_910", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_910", "test2 and state 911", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 911")

cmd.set_color("dummy_911", (0.8132256824298348, 0.9209534794309882, 0.9540945790080738))
cmd.color("dummy_911", "test2 and state 912", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 912")

cmd.set_color("dummy_912", (0.9983852364475202, 0.9499423298731258, 0.673125720876586))
cmd.color("dummy_912", "test2 and state 913", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 913")

cmd.set_color("dummy_913", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_913", "test2 and state 914", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 914")

cmd.set_color("dummy_914", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_914", "test2 and state 915", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 915")

cmd.set_color("dummy_915", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_915", "test2 and state 916", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 916")

cmd.set_color("dummy_916", (0.9977700884275279, 0.930872741253364, 0.6442137639369473))
cmd.color("dummy_916", "test2 and state 917", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 917")

cmd.set_color("dummy_917", (0.9167243367935409, 0.3430219146482122, 0.22399077277970011))
cmd.color("dummy_917", "test2 and state 918", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 918")

cmd.set_color("dummy_918", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_918", "test2 and state 919", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 919")

cmd.set_color("dummy_919", (0.9301038062283737, 0.3711649365628604, 0.23690888119953865))
cmd.color("dummy_919", "test2 and state 920", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 920")

cmd.set_color("dummy_920", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_920", "test2 and state 921", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 921")

cmd.set_color("dummy_921", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_921", "test2 and state 922", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 922")

cmd.set_color("dummy_922", (0.6240676662821992, 0.8137639369473281, 0.8934256055363322))
cmd.color("dummy_922", "test2 and state 923", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 923")

cmd.set_color("dummy_923", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_923", "test2 and state 924", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 924")

cmd.set_color("dummy_924", (0.43321799307958486, 0.6525951557093427, 0.8062283737024222))
cmd.color("dummy_924", "test2 and state 925", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 925")

cmd.set_color("dummy_925", (0.9991541714725106, 0.9737793156478277, 0.7092656670511341))
cmd.color("dummy_925", "test2 and state 926", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 926")

cmd.set_color("dummy_926", (0.9990003844675125, 0.9690119184928874, 0.7020376778162245))
cmd.color("dummy_926", "test2 and state 927", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 927")

cmd.set_color("dummy_927", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_927", "test2 and state 928", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 928")

cmd.set_color("dummy_928", (0.9993079584775086, 0.9785467128027683, 0.716493656286044))
cmd.color("dummy_928", "test2 and state 929", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 929")

cmd.set_color("dummy_929", (0.9928489042675894, 0.7169550173010381, 0.4129181084198385))
cmd.color("dummy_929", "test2 and state 930", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 930")

cmd.set_color("dummy_930", (0.831603229527105, 0.17716262975778546, 0.15271049596309114))
cmd.color("dummy_930", "test2 and state 931", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 931")

cmd.set_color("dummy_931", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_931", "test2 and state 932", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 932")

cmd.set_color("dummy_932", (0.5648596693579393, 0.7663975394079201, 0.8675893886966551))
cmd.color("dummy_932", "test2 and state 933", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 933")

cmd.set_color("dummy_933", (0.2742022299115725, 0.4631295655517109, 0.7081122645136487))
cmd.color("dummy_933", "test2 and state 934", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 934")

cmd.set_color("dummy_934", (0.9977700884275279, 0.930872741253364, 0.6442137639369473))
cmd.color("dummy_934", "test2 and state 935", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 935")

cmd.set_color("dummy_935", (0.9165705497885429, 0.9677047289504037, 0.9024221453287196))
cmd.color("dummy_935", "test2 and state 936", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 936")

cmd.set_color("dummy_936", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_936", "test2 and state 937", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 937")

cmd.set_color("dummy_937", (0.43321799307958486, 0.6525951557093427, 0.8062283737024222))
cmd.color("dummy_937", "test2 and state 938", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 938")

cmd.set_color("dummy_938", (0.9991541714725106, 0.9737793156478277, 0.7092656670511341))
cmd.color("dummy_938", "test2 and state 939", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 939")

cmd.set_color("dummy_939", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_939", "test2 and state 940", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 940")

cmd.set_color("dummy_940", (0.8944252210688197, 0.2961168781237985, 0.20246059207996925))
cmd.color("dummy_940", "test2 and state 941", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 941")

cmd.set_color("dummy_941", (0.615609381007305, 0.8069973087274126, 0.8897347174163783))
cmd.color("dummy_941", "test2 and state 942", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 942")

cmd.set_color("dummy_942", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_942", "test2 and state 943", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 943")

cmd.set_color("dummy_943", (0.6325259515570935, 0.8205305651672434, 0.8971164936562861))
cmd.color("dummy_943", "test2 and state 944", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 944")

cmd.set_color("dummy_944", (0.9957708573625529, 0.8630526720492119, 0.5502499038831219))
cmd.color("dummy_944", "test2 and state 945", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 945")

cmd.set_color("dummy_945", (0.9167243367935409, 0.3430219146482122, 0.22399077277970011))
cmd.color("dummy_945", "test2 and state 946", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 946")

cmd.set_color("dummy_946", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_946", "test2 and state 947", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 947")

cmd.set_color("dummy_947", (0.8213763936947329, 0.9249519415609382, 0.956401384083045))
cmd.color("dummy_947", "test2 and state 948", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 948")

cmd.set_color("dummy_948", (0.4404459823144945, 0.661207227989235, 0.8106881968473664))
cmd.color("dummy_948", "test2 and state 949", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 949")

cmd.set_color("dummy_949", (0.9070357554786621, 0.9640138408304498, 0.9199538638985004))
cmd.color("dummy_949", "test2 and state 950", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 950")

cmd.set_color("dummy_950", (0.6828143021914649, 0.8569780853517878, 0.9171856978085352))
cmd.color("dummy_950", "test2 and state 951", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 951")

cmd.set_color("dummy_951", (0.4887351018838909, 0.7054978854286814, 0.8343713956170704))
cmd.color("dummy_951", "test2 and state 952", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 952")

cmd.set_color("dummy_952", (0.9931564782775856, 0.7323337178008458, 0.42737408688965783))
cmd.color("dummy_952", "test2 and state 953", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 953")

cmd.set_color("dummy_953", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_953", "test2 and state 954", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 954")

cmd.set_color("dummy_954", (0.996078431372549, 0.8784313725490196, 0.5647058823529412))
cmd.color("dummy_954", "test2 and state 955", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 955")

cmd.set_color("dummy_955", (0.8975009611687813, 0.960322952710496, 0.9374855824682813))
cmd.color("dummy_955", "test2 and state 956", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 956")

cmd.set_color("dummy_956", (0.9122645136485967, 0.33364090734332946, 0.21968473663975396))
cmd.color("dummy_956", "test2 and state 957", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 957")

cmd.set_color("dummy_957", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_957", "test2 and state 958", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 958")

cmd.set_color("dummy_958", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_958", "test2 and state 959", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 959")

cmd.set_color("dummy_959", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_959", "test2 and state 960", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 960")

cmd.set_color("dummy_960", (0.9996155324875048, 0.988081507112649, 0.7309496347558632))
cmd.color("dummy_960", "test2 and state 961", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 961")

cmd.set_color("dummy_961", (0.9070357554786621, 0.9640138408304498, 0.9199538638985004))
cmd.color("dummy_961", "test2 and state 962", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 962")

cmd.set_color("dummy_962", (0.4802768166089966, 0.6987312572087659, 0.8306805074971165))
cmd.color("dummy_962", "test2 and state 963", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 963")

cmd.set_color("dummy_963", (0.9900807381776241, 0.6673587081891583, 0.3734717416378317))
cmd.color("dummy_963", "test2 and state 964", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 964")

cmd.set_color("dummy_964", (0.9070357554786621, 0.9640138408304498, 0.9199538638985004))
cmd.color("dummy_964", "test2 and state 965", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 965")

cmd.set_color("dummy_965", (0.5564013840830451, 0.7596309111880047, 0.8638985005767013))
cmd.color("dummy_965", "test2 and state 966", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 966")

cmd.set_color("dummy_966", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_966", "test2 and state 967", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 967")

cmd.set_color("dummy_967", (0.7724721261053442, 0.900961168781238, 0.942560553633218))
cmd.color("dummy_967", "test2 and state 968", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 968")

cmd.set_color("dummy_968", (0.31757016532103044, 0.5148019992310651, 0.7348712033833141))
cmd.color("dummy_968", "test2 and state 969", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 969")

cmd.set_color("dummy_969", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_969", "test2 and state 970", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 970")

cmd.set_color("dummy_970", (0.9942329873125721, 0.7861591695501731, 0.47797001153402535))
cmd.color("dummy_970", "test2 and state 971", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 971")

cmd.set_color("dummy_971", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_971", "test2 and state 972", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 972")

cmd.set_color("dummy_972", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_972", "test2 and state 973", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 973")

cmd.set_color("dummy_973", (0.9988465974625145, 0.9642445213379469, 0.6948096885813149))
cmd.color("dummy_973", "test2 and state 974", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 974")

cmd.set_color("dummy_974", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_974", "test2 and state 975", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 975")

cmd.set_color("dummy_975", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_975", "test2 and state 976", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 976")

cmd.set_color("dummy_976", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_976", "test2 and state 977", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 977")

cmd.set_color("dummy_977", (0.9118031526336026, 0.9658592848904267, 0.9111880046136099))
cmd.color("dummy_977", "test2 and state 978", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 978")

cmd.set_color("dummy_978", (0.9914648212226067, 0.6773548635140331, 0.3780853517877739))
cmd.color("dummy_978", "test2 and state 979", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 979")

cmd.set_color("dummy_979", (0.5986928104575164, 0.7934640522875818, 0.8823529411764706))
cmd.color("dummy_979", "test2 and state 980", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 980")

cmd.set_color("dummy_980", (0.9970011534025375, 0.9070357554786621, 0.6080738177623991))
cmd.color("dummy_980", "test2 and state 981", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 981")

cmd.set_color("dummy_981", (0.7398692810457518, 0.884967320261438, 0.9333333333333333))
cmd.color("dummy_981", "test2 and state 982", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 982")

cmd.set_color("dummy_982", (0.6240676662821992, 0.8137639369473281, 0.8934256055363322))
cmd.color("dummy_982", "test2 and state 983", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 983")

cmd.set_color("dummy_983", (0.9926951172625913, 0.7092656670511341, 0.4056901191849288))
cmd.color("dummy_983", "test2 and state 984", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 984")

cmd.set_color("dummy_984", (0.9946943483275663, 0.8092272202998847, 0.4996539792387543))
cmd.color("dummy_984", "test2 and state 985", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 985")

cmd.set_color("dummy_985", (0.9785467128027682, 0.9916955017301038, 0.7884659746251441))
cmd.color("dummy_985", "test2 and state 986", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 986")

cmd.set_color("dummy_986", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_986", "test2 and state 987", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 987")

cmd.set_color("dummy_987", (0.96239907727797, 0.46743560169165704, 0.281199538638985))
cmd.color("dummy_987", "test2 and state 988", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 988")

cmd.set_color("dummy_988", (0.9973087274125336, 0.9165705497885428, 0.6225297962322184))
cmd.color("dummy_988", "test2 and state 989", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 989")

cmd.set_color("dummy_989", (0.4971933871587852, 0.7122645136485968, 0.8380622837370243))
cmd.color("dummy_989", "test2 and state 990", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 990")

cmd.set_color("dummy_990", (0.40430603613994626, 0.6181468665897734, 0.7883890811226452))
cmd.color("dummy_990", "test2 and state 991", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 991")

cmd.set_color("dummy_991", (0.7969242599000386, 0.9129565551710881, 0.9494809688581315))
cmd.color("dummy_991", "test2 and state 992", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 992")

cmd.set_color("dummy_992", (0.9499423298731258, 0.9806228373702423, 0.8410611303344866))
cmd.color("dummy_992", "test2 and state 993", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 993")

cmd.set_color("dummy_993", (0.9951557093425606, 0.8322952710495963, 0.5213379469434832))
cmd.color("dummy_993", "test2 and state 994", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 994")

cmd.set_color("dummy_994", (0.7806228373702423, 0.904959630911188, 0.9448673587081892))
cmd.color("dummy_994", "test2 and state 995", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 995")

cmd.set_color("dummy_995", (0.44767397154940414, 0.6698193002691274, 0.8151480199923107))
cmd.color("dummy_995", "test2 and state 996", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 996")

cmd.set_color("dummy_996", (0.9356401384083045, 0.9750865051903114, 0.8673587081891581))
cmd.color("dummy_996", "test2 and state 997", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 997")

cmd.set_color("dummy_997", (0.964244521337947, 0.986159169550173, 0.8147635524798154))
cmd.color("dummy_997", "test2 and state 998", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 998")

cmd.set_color("dummy_998", (0.6663590926566706, 0.8475970780469051, 0.9118800461361015))
cmd.color("dummy_998", "test2 and state 999", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 999")

cmd.set_color("dummy_999", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_999", "test2 and state 1000", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1000")

cmd.set_color("dummy_1000", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_1000", "test2 and state 1001", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1001")

cmd.set_color("dummy_1001", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_1001", "test2 and state 1002", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1002")

cmd.set_color("dummy_1002", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_1002", "test2 and state 1003", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1003")

cmd.set_color("dummy_1003", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_1003", "test2 and state 1004", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1004")

cmd.set_color("dummy_1004", (0.9999231064975009, 0.9976163014225298, 0.7454056132256824))
cmd.color("dummy_1004", "test2 and state 1005", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1005")

cmd.set_color("dummy_1005", (0.7700884275278739, 0.11810841983852365, 0.15148019992310652))
cmd.color("dummy_1005", "test2 and state 1006", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1006")

cmd.set_color("dummy_1006", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_1006", "test2 and state 1007", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1007")

cmd.set_color("dummy_1007", (0.6409842368319878, 0.8272971933871589, 0.90080738177624))
cmd.color("dummy_1007", "test2 and state 1008", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1008")

cmd.set_color("dummy_1008", (0.972087658592849, 0.5374086889657824, 0.31349480968858134))
cmd.color("dummy_1008", "test2 and state 1009", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1009")

cmd.set_color("dummy_1009", (0.6991157247212612, 0.8649750096116878, 0.9217993079584775))
cmd.color("dummy_1009", "test2 and state 1010", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1010")

cmd.set_color("dummy_1010", (0.8702806612841217, 0.9489427143406383, 0.9702422145328721))
cmd.color("dummy_1010", "test2 and state 1011", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1011")

cmd.set_color("dummy_1011", (0.5141099577085737, 0.7257977700884276, 0.845444059976932))
cmd.color("dummy_1011", "test2 and state 1012", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1012")

cmd.set_color("dummy_1012", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_1012", "test2 and state 1013", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1013")

cmd.set_color("dummy_1013", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_1013", "test2 and state 1014", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1014")

cmd.set_color("dummy_1014", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_1014", "test2 and state 1015", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1015")

cmd.set_color("dummy_1015", (0.9737793156478277, 0.9898500576701268, 0.7972318339100347))
cmd.color("dummy_1015", "test2 and state 1016", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1016")

cmd.set_color("dummy_1016", (0.9930026912725874, 0.724644367550942, 0.4201460976547482))
cmd.color("dummy_1016", "test2 and state 1017", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1017")

cmd.set_color("dummy_1017", (0.6579008073817763, 0.8408304498269896, 0.9081891580161476))
cmd.color("dummy_1017", "test2 and state 1018", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1018")

cmd.set_color("dummy_1018", (0.5225682429834679, 0.7325643983083431, 0.8491349480968858))
cmd.color("dummy_1018", "test2 and state 1019", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1019")

cmd.set_color("dummy_1019", (0.9973087274125336, 0.9165705497885428, 0.6225297962322184))
cmd.color("dummy_1019", "test2 and state 1020", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1020")

cmd.set_color("dummy_1020", (0.988081507112649, 0.9953863898500577, 0.7709342560553633))
cmd.color("dummy_1020", "test2 and state 1021", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1021")

cmd.set_color("dummy_1021", (0.5733179546328336, 0.7731641676278356, 0.871280276816609))
cmd.color("dummy_1021", "test2 and state 1022", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1022")

cmd.set_color("dummy_1022", (0.8784313725490197, 0.9529411764705882, 0.9725490196078429))
cmd.color("dummy_1022", "test2 and state 1023", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1023")

cmd.set_color("dummy_1023", (0.8458285274894273, 0.9369473279507882, 0.9633217993079585))
cmd.color("dummy_1023", "test2 and state 1024", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1024")

cmd.set_color("dummy_1024", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_1024", "test2 and state 1025", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1025")

cmd.set_color("dummy_1025", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_1025", "test2 and state 1026", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1026")

cmd.set_color("dummy_1026", (0.9594771241830066, 0.9843137254901961, 0.8235294117647058))
cmd.color("dummy_1026", "test2 and state 1027", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1027")

cmd.set_color("dummy_1027", (0.6240676662821992, 0.8137639369473281, 0.8934256055363322))
cmd.color("dummy_1027", "test2 and state 1028", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1028")

cmd.set_color("dummy_1028", (0.5056516724336794, 0.7190311418685121, 0.8417531718569781))
cmd.color("dummy_1028", "test2 and state 1029", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1029")

cmd.set_color("dummy_1029", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_1029", "test2 and state 1030", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1030")

cmd.set_color("dummy_1030", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_1030", "test2 and state 1031", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1031")

cmd.set_color("dummy_1031", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_1031", "test2 and state 1032", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1032")

cmd.set_color("dummy_1032", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_1032", "test2 and state 1033", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1033")

cmd.set_color("dummy_1033", (0.4404459823144945, 0.661207227989235, 0.8106881968473664))
cmd.color("dummy_1033", "test2 and state 1034", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1034")

cmd.set_color("dummy_1034", (0.9999231064975009, 0.9976163014225298, 0.7454056132256824))
cmd.color("dummy_1034", "test2 and state 1035", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1035")

cmd.set_color("dummy_1035", (0.9985390234525182, 0.9547097270280661, 0.6803537101114956))
cmd.color("dummy_1035", "test2 and state 1036", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1036")

cmd.set_color("dummy_1036", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_1036", "test2 and state 1037", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1037")

cmd.set_color("dummy_1037", (0.7154171472510573, 0.8729719338715878, 0.9264129181084199))
cmd.color("dummy_1037", "test2 and state 1038", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1038")

cmd.set_color("dummy_1038", (0.8879661668589005, 0.9566320645905421, 0.955017301038062))
cmd.color("dummy_1038", "test2 and state 1039", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1039")

cmd.set_color("dummy_1039", (0.9979238754325259, 0.9356401384083045, 0.6514417531718569))
cmd.color("dummy_1039", "test2 and state 1040", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1040")

cmd.set_color("dummy_1040", (0.9261053440984237, 0.9713956170703576, 0.8848904267589387))
cmd.color("dummy_1040", "test2 and state 1041", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1041")

cmd.set_color("dummy_1041", (0.993925413302576, 0.7707804690503652, 0.463514033064206))
cmd.color("dummy_1041", "test2 and state 1042", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1042")

cmd.set_color("dummy_1042", (0.8621299500192235, 0.9449442522106882, 0.9679354094579009))
cmd.color("dummy_1042", "test2 and state 1043", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1043")

cmd.set_color("dummy_1043", (0.38262206843521723, 0.5923106497500962, 0.7750096116878125))
cmd.color("dummy_1043", "test2 and state 1044", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1044")

cmd.set_color("dummy_1044", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_1044", "test2 and state 1045", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1045")

cmd.set_color("dummy_1045", (0.5902345251826222, 0.7866974240676664, 0.8786620530565168))
cmd.color("dummy_1045", "test2 and state 1046", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1046")

cmd.set_color("dummy_1046", (0.9167243367935409, 0.3430219146482122, 0.22399077277970011))
cmd.color("dummy_1046", "test2 and state 1047", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1047")

cmd.set_color("dummy_1047", (0.8239138792772011, 0.16978085351787775, 0.15255670895809306))
cmd.color("dummy_1047", "test2 and state 1048", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1048")

cmd.set_color("dummy_1048", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_1048", "test2 and state 1049", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1049")

cmd.set_color("dummy_1049", (0.9943867743175702, 0.7938485198000771, 0.4851980007689352))
cmd.color("dummy_1049", "test2 and state 1050", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1050")

cmd.set_color("dummy_1050", (0.8831987697039602, 0.9547866205305652, 0.9637831603229524))
cmd.color("dummy_1050", "test2 and state 1051", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1051")

cmd.set_color("dummy_1051", (0.8831987697039602, 0.9547866205305652, 0.9637831603229524))
cmd.color("dummy_1051", "test2 and state 1052", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1052")

cmd.set_color("dummy_1052", (0.9651672433679355, 0.48742791234140714, 0.29042675893886966))
cmd.color("dummy_1052", "test2 and state 1053", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1053")

cmd.set_color("dummy_1053", (0.6071510957324107, 0.8002306805074972, 0.8860438292964244))
cmd.color("dummy_1053", "test2 and state 1054", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1054")

cmd.set_color("dummy_1054", (0.8050749711649368, 0.9169550173010381, 0.9517877739331027))
cmd.color("dummy_1054", "test2 and state 1055", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1055")

cmd.set_color("dummy_1055", (0.9991541714725106, 0.9737793156478277, 0.7092656670511341))
cmd.color("dummy_1055", "test2 and state 1056", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1056")

cmd.set_color("dummy_1056", (0.9499423298731258, 0.9806228373702423, 0.8410611303344866))
cmd.color("dummy_1056", "test2 and state 1057", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1057")

cmd.set_color("dummy_1057", (0.9963860053825452, 0.8879661668589004, 0.5791618608227604))
cmd.color("dummy_1057", "test2 and state 1058", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1058")

cmd.set_color("dummy_1058", (0.9859284890426759, 0.6373702422145329, 0.3596309111880046))
cmd.color("dummy_1058", "test2 and state 1059", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1059")

cmd.set_color("dummy_1059", (0.690965013456363, 0.8609765474817378, 0.9194925028835064))
cmd.color("dummy_1059", "test2 and state 1060", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1060")

cmd.set_color("dummy_1060", (0.6701268742791234, 0.02214532871972319, 0.14948096885813147))
cmd.color("dummy_1060", "test2 and state 1061", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1061")

cmd.set_color("dummy_1061", (0.7317185697808536, 0.880968858131488, 0.9310265282583622))
cmd.color("dummy_1061", "test2 and state 1062", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1062")

cmd.set_color("dummy_1062", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_1062", "test2 and state 1063", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1063")

cmd.set_color("dummy_1063", (0.9213379469434834, 0.9695501730103806, 0.8936562860438291))
cmd.color("dummy_1063", "test2 and state 1064", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1064")

cmd.set_color("dummy_1064", (0.9930026912725874, 0.724644367550942, 0.4201460976547482))
cmd.color("dummy_1064", "test2 and state 1065", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1065")

cmd.set_color("dummy_1065", (0.6746635909265668, 0.8529796232218377, 0.914878892733564))
cmd.color("dummy_1065", "test2 and state 1066", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1066")

cmd.set_color("dummy_1066", (0.7480199923106499, 0.888965782391388, 0.9356401384083045))
cmd.color("dummy_1066", "test2 and state 1067", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1067")

cmd.set_color("dummy_1067", (0.892733564013841, 0.9584775086505191, 0.9462514417531717))
cmd.color("dummy_1067", "test2 and state 1068", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1068")

cmd.set_color("dummy_1068", (0.9925413302575933, 0.7015763168012303, 0.3984621299500192))
cmd.color("dummy_1068", "test2 and state 1069", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1069")

cmd.set_color("dummy_1069", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_1069", "test2 and state 1070", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1070")

cmd.set_color("dummy_1070", (0.6494425221068819, 0.8340638216070742, 0.9044982698961938))
cmd.color("dummy_1070", "test2 and state 1071", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1071")

cmd.set_color("dummy_1071", (0.5394848135332565, 0.7460976547481738, 0.8565167243367935))
cmd.color("dummy_1071", "test2 and state 1072", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1072")

cmd.set_color("dummy_1072", (0.7887735486351405, 0.9089580930411381, 0.9471741637831603))
cmd.color("dummy_1072", "test2 and state 1073", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1073")

cmd.set_color("dummy_1073", (0.9803921568627452, 0.5973856209150327, 0.3411764705882353))
cmd.color("dummy_1073", "test2 and state 1074", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1074")

cmd.set_color("dummy_1074", (0.5648596693579393, 0.7663975394079201, 0.8675893886966551))
cmd.color("dummy_1074", "test2 and state 1075", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1075")

cmd.set_color("dummy_1075", (0.700884275278739, 0.05167243367935409, 0.1500961168781238))
cmd.color("dummy_1075", "test2 and state 1076", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1076")

cmd.set_color("dummy_1076", (0.2742022299115725, 0.4631295655517109, 0.7081122645136487))
cmd.color("dummy_1076", "test2 and state 1077", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1077")

cmd.set_color("dummy_1077", (0.9997693194925029, 0.9928489042675894, 0.7381776239907728))
cmd.color("dummy_1077", "test2 and state 1078", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1078")

cmd.set_color("dummy_1078", (0.9977700884275279, 0.930872741253364, 0.6442137639369473))
cmd.color("dummy_1078", "test2 and state 1079", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1079")

cmd.set_color("dummy_1079", (0.9951557093425606, 0.8322952710495963, 0.5213379469434832))
cmd.color("dummy_1079", "test2 and state 1080", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1080")

cmd.set_color("dummy_1080", (0.8295271049596311, 0.9289504036908882, 0.9587081891580161))
cmd.color("dummy_1080", "test2 and state 1081", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1081")

cmd.set_color("dummy_1081", (0.9707035755478662, 0.5274125336409073, 0.308881199538639))
cmd.color("dummy_1081", "test2 and state 1082", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1082")

cmd.set_color("dummy_1082", (0.9690119184928874, 0.9880046136101499, 0.805997693194925))
cmd.color("dummy_1082", "test2 and state 1083", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1083")

cmd.set_color("dummy_1083", (0.9022683583237218, 0.962168396770473, 0.9287197231833908))
cmd.color("dummy_1083", "test2 and state 1084", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1084")

cmd.set_color("dummy_1084", (0.9845444059976932, 0.6273740868896578, 0.3550173010380623))
cmd.color("dummy_1084", "test2 and state 1085", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1085")

cmd.set_color("dummy_1085", (0.5817762399077278, 0.7799307958477508, 0.8749711649365628))
cmd.color("dummy_1085", "test2 and state 1086", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1086")

cmd.set_color("dummy_1086", (0.8376778162245292, 0.9329488658208382, 0.9610149942329873))
cmd.color("dummy_1086", "test2 and state 1087", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1087")

cmd.set_color("dummy_1087", (0.3898500576701269, 0.6009227220299886, 0.7794694348327567))
cmd.color("dummy_1087", "test2 and state 1088", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1088")

cmd.set_color("dummy_1088", (0.7235678585159555, 0.8769703960015379, 0.928719723183391))
cmd.color("dummy_1088", "test2 and state 1089", 1)

cmd.set("stick_transparency", 0.5, f"test2 and state 1089")

for x in positions_viewport_callbacks:
    for y in positions_viewport_callbacks[x]:
        positions_viewport_callbacks[x][y].load()

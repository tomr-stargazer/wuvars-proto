""" 
This is a helper file to `cygob7_jjh_script.py`.

It defines things that I'd rather not have stealing space in
the script proper, because they take up tons of space.

"""

import numpy as np

wise_disks = np.array([44027709956392,
44027709959743,
44027710034711,
-1,
44027710079437,
44027709895213,
44027709997118,
44027709993435,
44027710056350,
44027709899303,
44027710021056,
44027710097045,
44027710064752,
44027710026304])

wise_disks_names = np.array([14537,
16393,
22088,
29302,
31002,
31750,
31848,
32627,
33466,
34643,
34726,
35166,
35254,
35280])

wise_extras = np.array([44027709957958,
44027709999438,
44027710155529,
-1,
44027709895321,
44027709974258,
44027709894711,
44027710027179,
44027710053494,
44027710051061,
44027710018938,
44027709894777,
44027710021888,
44027709985940,
44027710079486])

wise_extras_names = np.array([12017,
15170,
27251,
29302,
29462,
29815,
30034,
30943,
32121,
33383,
33479,
33909,
34276,
34655,
34904])

aspin_sources = np.array([-1,
44027709895323,
44027709926268,
44027709926317,
44027710047394,
44027709871330,
44027709925375,
44027709986782,
44027710046074,
44027710046077])

aspin_sources_names = np.array(['aspin Braid star',
'aspin CN 1S',
'aspin CN 3N',
'aspin CN 6N',
'aspin Cyg 19',
'aspin HH381 IRS ',
'aspin HH627-STAR',
'aspin IRAS 14',
'aspin IRAS 15N',
'aspin IRAS 15S'])

wise_trans = np.array([44027709870583,
44027710026767,
44027710026681,
44027710025882,
44027710025862,
44027710025741,
44027709872688,
44027709872058,
44027709871821,
44027710022475,
44027710022381,
-1,
44027709869803,
-1,
44027710020232,
44027710020128,
44027710057220,
-1,
44027710054362,
44027710053092,
44027710052979,
44027710053089,
-1,
-1,
-1,
-1,
-1,
44027709925375,
-1,
-1,
-1,
44027710046920,
44027710089066,
44027709895323,
44027710080715,
44027710081177,
44027709978489,
44027709984925,
44027709784533,
-1,
-1,
-1,
44027709978531,
44027709907407,
44027709840086,
44027709998214,
44027709998291,
44027710110481,
-1,
-1,
44027710074434,
44027710074374,
44027710113810,
-1,
-1,
-1,
44027710004926,
44027710000081,
44027710006032,
44027709958000])

wise_trans_names = np.array([23626,
31690,
30603,
31170,
31666,
33482,
22443,
34493,
32237,
34473,
33711,
30607,
29897,
30942,
35257,
31673,
32723,
35252,
30918,
32629,
33631,
31586,
31643,
7377,
5212,
8845,
29024,
34305,
7198,
9295,
8297,
32637,
28616,
34808,
27320,
28614,
29025,
30173,
15184,
7698,
9304,
6505,
30833,
12775,
13997,
12266,
13323,
26269,
3414,
8574,
32340,
28997,
28576,
7736,
8562,
8025,
10951,
11810,
16077,
15289])

welo_sources = np.array([44027709894711, 44027709894777])

welo_names = wise_extras_names[ np.in1d(wise_extras, welo_sources) ]

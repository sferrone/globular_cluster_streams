import json
from create_interval_json import create_json

# 14 may
# before 3 may
create_json("svol",time=[-1,1],lon=[34,53],lat=[26,85],D=[0,9],PMB=[-5,0],PML_COSB=[-5,0])
create_json("gaia-1",time=[-1,1],lon=[-68,-29],lat=[40,70],D=[1,5],PMB=[-25,-6],PML_COSB=[-25,-7])
create_json("gaia-2",time=[-1,1],lon=[87,97],lat=[-88,-70],D=[2,12],PMB=[-7,4],PML_COSB=[-7,4])
create_json("indus",time=[-1,1],lon=[-60,-10],lat=[-60,-40],D=[10,30],PMB=[0,5],PML_COSB=[-10,-5])
create_json("m5",time=[-1,1],lon=[-55,10],lat=[40,85],D=[3,12],PMB=[-15,-5],PML_COSB=[-8,6])
create_json("leiptr",time=[-1,1],lon=[-178,-92],lat=[-45,-5],D=[5,15],PMB=[0,9],PML_COSB=[9,18])
create_json("NGC6397",time=[-1,1],lon=[-28,-10],lat=[-30,-5],D=[0,7],PMB=[-30,30],PML_COSB=[-30,30])
# updates from earding the paper
create_json("gaia-7",time=[-1,1],lon=[-90,-60],lat=[40,55],D=[3, 8],PMB=[-6,1],PML_COSB=[-12,-5])
create_json("gaia-10",time=[-1,1],lon=[-145,-115],lat=[45,65],D=[10,20],PMB=[-9,-4],PML_COSB=[1,5])


# I think the proper motions are in RA and dec for some of the paper yet in L and B for some other parts of the paper




create_json("fjorm",time=[-1,1],lon=[ -70, 110], lat=[30, 89], D=[1, 5], PMB=[-8, 8], PML_COSB=[-5, 11])
create_json("gjoll",    time=[-1,1],lon=[-177,-110], lat=[-34,-4],  D=[0, 9],  PMB=[  3, 19], PML_COSB=[15, 30])
create_json("fimbulthul",time=[-1,1],lon=[-25,-63],lat=[20,40],D=[0,8],PMB=[-13,-2],PML_COSB=[-30,0])
create_json("NGC288",time=[-1,1],lon=[-180, 180], lat=[-90,-80], D=[0,13], PMB=[-10, 10],  PML_COSB=[-2, 10])
create_json("NGC1261",time=[-1,1],lon=[-110, -70], lat=[-60,-30], D=[10,30], PMB=[ 0,  5],  PML_COSB=[-5, 5])
create_json("pal5",         time=[-1,1],lon=[ -15,  15], lat=[ 30, 50], D=[15, 32], PMB=[ 4,  5], PML_COSB=[-2, -7])
create_json("NGC1851",time=[-1,1],lon=[-130,-100], lat=[-45,-20], D=[10,30], PMB=[ 1,  5],  PML_COSB=[-1, 3])



create_json("slidr",time=[-1,1],lon=[-150, -85], lat=[38, 68], D=[1, 5], PMB=[-25, -8], PML_COSB=[ -15, -8])
create_json("m92",  time=[-1,1],lon=[  61,  75], lat=[ 26, 40], D=[ 5,13], PMB=[  3,  7], PML_COSB=[-3, 2])
create_json("orphan",time=[-1,1],lon=[-60,-180],lat=[20,55],D=[10,30],PMB=[-5,5],PML_COSB=[-10,2])
create_json("phlegethon",time=[-1,1],lon=[-15,80],lat=[-25,-55],D=[0,5],PMB=[-19,10],PML_COSB=[-15,-30])
create_json("NGC3201",      time=[-1,1],lon=[ -89, -75], lat=[  5, 15], D=[ 1, 5], PMB=[  2, 8], PML_COSB=[-30, 30])
create_json("c-9a",     time=[-1,1],lon=[ 125, 180], lat=[-88,-70], D=[ 5,11], PMB=[-30,30], PML_COSB=[-30, 30])
create_json("c-9b",     time=[-1,1],lon=[-180,-130], lat=[-88,-70], D=[ 5,11], PMB=[-30,30], PML_COSB=[-30, 30])

create_json("phoenix",time=[-1,1],lon=[-110,-60],lat=[-75,-55],D=[15,20],PMB=[0,4],PML_COSB=[-5,0]) # balbinot et al
create_json("acs",time=[-1,1],lon=[-175,-120],lat=[20,40], D=[3, 12], PMB=[-2,4], PML_COSB=[-5,2])
create_json("NGC5466",time=[-1,1],lon=[22,50],lat=[52,83],D=[10,30],PMB=[5,15],PML_COSB=[-5,5])
create_json("c-11",time=[-1,1],lon=[15,45],lat=[-55,-35], D=[3, 12], PMB=[-3,3], PML_COSB=[-10,-3])


create_json("c-1",  time=[-1,1],lon=[ -50, -10], lat=[70, 80], D=[4, 9], PMB=[-15, -5], PML_COSB=[-7, -3])
create_json("c-12",time=[-1,1],lon=[90, 120],lat=[-55,-15],D=[10,30], PMB=[-5,-1],  PML_COSB=[-1,4])
create_json("c-5",          time=[-1,1],lon=[-180,-165],lat=[20,35],D=[1,5],PMB=[-7,4],PML_COSB=[6,30])
create_json("c-7", time=[-1,1],lon=[0,-30], lat=[-35,-5], D=[ 0, 6], PMB=[ 4,12], PML_COSB=[-13, -20])
create_json("c-8", time=[-1,1],lon=[-60,-35],lat=[-55,-31],D=[2,5],PMB=[-12,-7],PML_COSB=[-7,1])
create_json("c-10",time=[-1,1],lon=[155,175],lat=[17,42],D=[0,12],PMB=[5,11], PML_COSB=[1,8])
create_json("gaia-6",time=[-1,1],lon=[25,69],lat=[55,82],D=[5,12],PMB=[7,20],PML_COSB=[-10,1])
create_json("gaia-12",time=[-1,1],lon=[150,170],lat=[-50,-30],D=[8,17],PMB=[-3,3],PML_COSB=[8,15])

create_json("c-3",          time=[-1,1],lon=[ -37, -75], lat=[ 15, 45], D=[ 4,  8], PMB=[-11, -6], PML_COSB=[-14, -7])
create_json("gaia-8",time=[-1,1],lon=[-85,-40],lat=[20,55],D=[3,12],PMB=[-12,-5],PML_COSB=[-12,-2])

create_json("c-2",  time=[-1,1],lon=[70, 120], lat=[ 40, 50], D=[ 1, 6], PMB=[1, 8], PML_COSB=[1, 8])
create_json("gaia-9",time=[-1,1],lon=[55,125],lat=[35,55],D=[2,8],PMB=[1,10],PML_COSB=[3,20])

create_json("c-6", time=[-1,1],lon=[50,85], lat=[ 7, 31], D=[ 7,15], PMB=[2, 9], PML_COSB=[2, 9])
create_json("gaia-11",time=[-1,1],lon=[60,90],lat=[15,40],D=[7,15],PMB=[5,10],PML_COSB=[2,8])

create_json("GD1a",time=[-1,1],lon=[180,90],lat=[26,70],D=[4,10],PMB=[-5,12],PML_COSB=[5,28]) # # This one is special because it crosses 180
create_json("GD1b",time=[-1,1],lon=[-180,-130],lat=[26,70],D=[4,10],PMB=[-10,-5],PML_COSB=[8,28]) 

create_json("gunnthra",time=[-1,1],lon=[-50,-15],lat=[-50,-20],D=[0,6],PMB=[7,25],PML_COSB=[-28,-15])
create_json("hrid",time=[-1,1],lon=[40,120],lat=[0,30],D=[0,7],PMB=[7,30],PML_COSB=[11,30]) # there are proper motions reported in the paper
create_json("jhelum",time=[-1,1],lon=[-5,-60],lat=[-40,-75],D=[5,13],PMB=[-9,2],PML_COSB=[-6,-14])
create_json("kshir",time=[-1,1],lon=[89,130],lat=[35,55],D=[7,13],PMB=[2,10],PML_COSB=[2,10])
create_json("sylgr",time=[-1,1],lon=[ -110, -60], lat=[42, 75], D=[0, 7], PMB=[-40, -8], PML_COSB=[ -6,-13])
create_json("vid",time=[-1,1],lon=[-180,-120],lat=[-90,-65],D=[15,30],PMB=[-8,3],PML_COSB=[-8,3])
create_json("ylgr",time=[-1,1],lon=[-124,-60],lat=[10,60],D=[5,16],PMB=[-3,-12], PML_COSB=[-1, 9])
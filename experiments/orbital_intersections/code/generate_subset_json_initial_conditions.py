import json
from create_interval_json import create_json
## FIX THESE
create_json("c-1",  time=[-1,1],lon=[ -50, -10], lat=[70, 80], D=[4, 9], PMB=[-15, -5], PML_COSB=[-7, -3])
create_json("c-12",time=[-1,1],lon=[90, 120],lat=[-55,-15],D=[10,30], PMB=[-5,-1],  PML_COSB=[-1,4])
create_json("c-2",  time=[-1,1],lon=[70, 120], lat=[ 40, 50], D=[ 1, 6], PMB=[1, 8], PML_COSB=[1, 8])
create_json("c-3",          time=[-1,1],lon=[ -37, -75], lat=[ 15, 45], D=[ 4,  8], PMB=[-11, -6], PML_COSB=[-14, -7])
create_json("c-5",          time=[-1,1],lon=[-180,-165],lat=[20,35],D=[1,5],PMB=[-7,4],PML_COSB=[6,30])
create_json("c-6",  time=[-1,1],lon=[50,85], lat=[ 7, 31], D=[ 7,15], PMB=[2, 9], PML_COSB=[2, 9])
create_json("c-7", time=[-1,1],lon=[0,-30], lat=[-35,-5], D=[ 0, 6], PMB=[ 4,12], PML_COSB=[-13, -20])
create_json("c-8", time=[-1,1],lon=[-60,-35],lat=[-55,-31],D=[2,5],PMB=[-12,-7],PML_COSB=[-7,1])
create_json("c-9a",     time=[-1,1],lon=[ 125, 180], lat=[-88,-70], D=[ 0,11], PMB=[-7,1], PML_COSB=[-7,1])
create_json("c-9b",     time=[-1,1],lon=[-180,-130], lat=[-88,-70], D=[ 0,11], PMB=[-7,1], PML_COSB=[-7,1])
create_json("c-10",time=[-1,1],lon=[155,175],lat=[17,42],D=[0,12],PMB=[5,11], PML_COSB=[1,8])
create_json("gaia-6",time=[-1,1],lon=[25,69],lat=[55,82],D=[5,12],PMB=[7,20],PML_COSB=[-10,1])
create_json("gaia-9",time=[-1,1],lon=[55,125],lat=[35,55],D=[2,8],PMB=[1,10],PML_COSB=[3,20])
create_json("GD1a",time=[-1,1],lon=[180,90],lat=[26,70],D=[4,10],PMB=[-5,12],PML_COSB=[5,28]) # # This one is special because it crosses 180
create_json("GD1b",time=[-1,1],lon=[-180,-130],lat=[26,70],D=[4,10],PMB=[-10,-5],PML_COSB=[8,28]) 
create_json("gunnthra",time=[-1,1],lon=[-50,-15],lat=[-50,-20],D=[0,6],PMB=[7,25],PML_COSB=[-28,-15])
create_json("hrid",time=[-1,1],lon=[40,120],lat=[0,30],D=[0,7],PMB=[7,30],PML_COSB=[11,30]) # there are proper motions reported in the paper
create_json("jhelum",time=[-1,1],lon=[-5,-60],lat=[-40,-75],D=[5,13],PMB=[-9,2],PML_COSB=[-6,-14])
create_json("kshir",time=[-1,1],lon=[89,130],lat=[35,55],D=[7,13],PMB=[2,10],PML_COSB=[2,10])
create_json("sylgr",time=[-1,1],lon=[-110,-60],lat=[42,75],D=[0,7],PMB=[-40,-8],PML_COSB=[-6,-13])
create_json("vid",time=[-1,1],lon=[-180,-120],lat=[-90,-65],D=[15,30],PMB=[-8,3],PML_COSB=[-8,3])
create_json("ylgr",time=[-1,1],lon=[-124,-60],lat=[10,60],D=[5,16],PMB=[-3,-12], PML_COSB=[-1, 9])
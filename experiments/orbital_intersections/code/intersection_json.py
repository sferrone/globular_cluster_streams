import h5py
import json
import numpy as np 
import astropy.units as u
import astropy.coordinates as coords
from astropy.coordinates import Galactic
import sys 
sys.path.append('../../functions')
import global_functions as my_gf 

# Requisites are list of clusters
cluster_fname   = "mini_list_of_clusters.txt"
# path to json object with the intervalz
interval_fname  = sys.argv[1]
# set path information to the orbits of each cluster
ftype           = "backwardorbits"
back_path       = "../the_orbits/backward"
forward_path    = "../the_orbits/forward"
output_path     = "../outputs/" 

# get intersection of desired clusters and those in the avaliable in the directory
clusters, n_yeras   = my_gf.match_cluster_list_to_data_directory(cluster_fname, forward_path, ftype)
fp                  = open(output_path+interval_fname+".json", "r+")
interval            = json.load(fp)

# Transform JSON object into appropriate units    
interval['LONG']        = interval['LONG']      * u.degree
interval['LAT']         = interval['LAT']       * u.degree
interval['D']           = interval['D']         * u.kpc
interval['PMB']         = interval['PMB']       * u.mas/u.year
interval['PML_COSB']    = interval['PML_COSB']  * u.mas/u.year
# what about the line of sight velocity ?
# There's no radial velocity in streams under 'galactic'. But equatorial has both RA and RV? 
# initialize output of matches
matches_dict = {}
n_clusters = len(clusters)
print("Doing {:d} clusters for {:s}".format(int(n_clusters), interval_fname))
for i in range(len(clusters)):
    progress = np.floor(100*i/n_clusters)
    if (progress % 10 == 0):
        print("{:d} %".format(int(progress)))    
    # open forward and backward orbits
    tb,xb,yb,zb,vxb,vyb,vzb = np.loadtxt(back_path+"/"+"orbit"+clusters[i]+".dat", unpack=True,skiprows=1)
    tf,xf,yf,zf,vxf,vyf,vzf = np.loadtxt(forward_path+"/"+"orbit"+clusters[i]+".dat", unpack=True,skiprows=1)
    # TIME was given as 5 GYRs ago <- 0. So it starts a 5 ends at today
    tb = -np.flip(tb)  # Now in billions of years and moves forward
    xb = np.flip(xb)  
    yb = np.flip(yb)  
    zb = np.flip(zb)  
    vxb = -np.flip(vxb)
    vyb = -np.flip(vyb)
    vzb = -np.flip(vzb) 

    # concatenate backward and forward 
    t = np.hstack([tb,tf])/10
    x = np.hstack([xb,xf]) * u.kpc
    y = np.hstack([yb,yf]) * u.kpc
    z = np.hstack([zb,zf]) * u.kpc
    # Units were given in 10 km / s (hectoras anesantian units). Put them in km / s
    vx = np.hstack([vxb,vxf]) * u.km / u.s * 10
    vy = np.hstack([vyb,vyf]) * u.km / u.s * 10
    vz = np.hstack([vzb,vzf])  * u.km / u.s * 10
    
    # put object into skycoord object
    galacto_centric = coords.Galactocentric(x=x,y=y,z=z,v_x=vx,v_y=vy, v_z=vz)
    # transform to galactic coordinates
    galactic_coordinates = galacto_centric.transform_to(Galactic)
    my_lon_coords = galactic_coordinates.l.wrap_at(180*u.degree)
    # do the positional the criteria
    lon_criteria        = (interval['LONG'][0]  < my_lon_coords)                  & (my_lon_coords                  < interval['LONG'][1])
    lat_criteria        = (interval['LAT'][0]   < galactic_coordinates.b)         & (galactic_coordinates.b         < interval['LAT'][1])
    D_critera           = (interval["D"][0]     < galactic_coordinates.distance)  & (galactic_coordinates.distance  < interval["D"][1])
    positional_criteria = lon_criteria & lat_criteria & D_critera
    # do the velocity criteria
    PMB_criteria        = (interval["PMB"][0]       < galactic_coordinates.pm_b)      & (galactic_coordinates.pm_b        < interval["PMB"][1])
    PML_COSB_criteria   = (interval["PML_COSB"][0]  < galactic_coordinates.pm_l_cosb) & (galactic_coordinates.pm_l_cosb   < interval["PML_COSB"][1])
    velocity_criteria   = PMB_criteria & PML_COSB_criteria
    # add time limit 
    time_critera        = (interval["time"][0]  < t) & (t   < interval["time"][1])
    # final criterian
    criterian           = time_critera & velocity_criteria & positional_criteria
    # store result into python object only if there are stars
    if np.sum(criterian) > 0:
        matches_dict[clusters[i]] = t[criterian].tolist()
# strip of units to prepare for json output
for key in interval.keys():
    if (key != "matches") and (key != "time") :
        interval[key] = interval[key].value
        interval[key] = interval[key].tolist()
if bool(matches_dict):
    interval['matches'] = matches_dict
else:   
    interval["matches"] = "None"
# rewrite the interval file
fp.seek(0)
json.dump(interval, fp)
fp.close()
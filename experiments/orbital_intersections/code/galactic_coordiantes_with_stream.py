import h5py
import json
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import astropy.units as u
import astropy.coordinates as coords
from astropy.coordinates import Galactic
from astropy.io import fits
import sys 
sys.path.append('../../functions')
import global_functions as my_gf 

# Requisites are list of clusters
cluster  = sys.argv[1]
# path to json object with the intervalz
stream_name     = sys.argv[2]
interval_fname  = "../outputs/"+stream_name+".json"
# set path information to the orbits of each cluster
ftype           = "backwardorbits"
back_path       = "../the_orbits/backward"
forward_path    = "../the_orbits/forward"
output_path     = "../outputs/" 
# get intersection of desired clusters and those in the avaliable in the directory
fp                  = open(interval_fname, "r+")
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

# open forward and backward orbits
tb,xb,yb,zb,vxb,vyb,vzb = np.loadtxt(back_path+"/"+"orbit"+cluster+".dat", unpack=True,skiprows=1)
tf,xf,yf,zf,vxf,vyf,vzf = np.loadtxt(forward_path+"/"+"orbit"+cluster+".dat", unpack=True,skiprows=1)
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
vx = np.hstack([vxb,vxf]) * 10 * u.km / u.s 
vy = np.hstack([vyb,vyf]) * 10 * u.km / u.s 
vz = np.hstack([vzb,vzf]) * 10 * u.km / u.s

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
net_criterian       = time_critera & velocity_criteria & positional_criteria

# get the position of the mass today
index_cluster_today = np.argmin(abs(t))

# do plot of each quantity 
q_names = ["PML_COSB", "time", "PMB", "distance"]
for quantity_name in q_names:
    criterian_name = "time"
    dot_size = .5
    my_cmap = "jet_r"

    # set criteria for each 
    if criterian_name =="time":
        criterian = time_critera
    if quantity_name=="PML_COSB":
        vmin = -20
        vmax =  20
        my_quantity         = galactic_coordinates.pm_l_cosb[criterian]
        good_quantity       = galactic_coordinates.pm_l_cosb[net_criterian]
    elif quantity_name =="PMB":
        vmin = -20
        vmax =  20
        my_quantity         = galactic_coordinates.pm_b[criterian]
        good_quantity       = galactic_coordinates.pm_b[net_criterian]
    elif quantity_name =="time":
        my_quantity         = t[criterian]
        good_quantity       = t[net_criterian]
        vmin = -1
        vmax =  1
    elif quantity_name =="distance":
        my_quantity         = galactic_coordinates.distance[criterian]
        good_quantity       = galactic_coordinates.distance[net_criterian]
        vmin = 1.5
        vmax = 20
        my_cmap = "viridis"
    # get coordinates
    cbar_name           = quantity_name
    my_x_coordinate     = my_lon_coords[criterian]
    my_y_coordinate     = galactic_coordinates.b[criterian]
    good_x_coordinates  = my_lon_coords[net_criterian]
    good_y_coordinates  = galactic_coordinates.b[net_criterian]

    # start plot 
    my_fontsize = 20
    fig,ax = plt.subplots(1,1)
    fig.set_figheight(5)
    fig.set_figwidth(11)
    # create the shape showing the latitude and longitude range
    lat_width = interval['LAT'][1] - interval['LAT'][0]
    lon_width = interval['LONG'][1] - interval['LONG'][0]
    me_patch=patches.Rectangle((interval['LONG'][0].value, interval['LAT'][0].value), lon_width.value, lat_width.value, alpha=.4)
    im = ax.scatter(my_x_coordinate, my_y_coordinate, c=my_quantity, s=dot_size, zorder=2,cmap=my_cmap, vmin=vmin, vmax=vmax)
    # make the dots bigger if the criteria was satisfied
    if np.sum(net_criterian) > 0:
        im_good = ax.scatter(good_x_coordinates,good_y_coordinates,c=good_quantity,s=20*dot_size,zorder=3,cmap=my_cmap, vmin=vmin, vmax=vmax)
    # plt the cluster from today
    plt.scatter(my_lon_coords[index_cluster_today],galactic_coordinates.b[index_cluster_today], marker='o', c='k', zorder=4, alpha=.5)
    # make plot nice
    ax.add_patch(me_patch)
    ax.set_xlim(-180,180)
    ax.set_ylim(-90,90)
    ax.set_xticks(np.arange(-180,210,30))
    ax.set_yticks(np.arange(-90,120,30))
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_aspect('equal')
    cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
    cbar.ax.get_yaxis().labelpad = 20
    # add information 
    ax.set_title(stream_name+"-"+cluster, size=my_fontsize)
    ax.set_xlabel("LONG", size=my_fontsize)
    ax.set_ylabel("LAT",size=my_fontsize)
    ax.tick_params(labelsize=.75*my_fontsize)
    xpos = 180 - .05*(360)
    ypos = -90 + .02*(180)
    words = "PMB: [{:.1f} {:.1f}] mas/yr\nPML: [{:.1f} {:.1f}] mas/yr\nD:   [{:.1f} {:.1f}] kpc".format(interval['PMB'][0].value,interval["PMB"][1].value, interval["PML_COSB"][0].value, interval["PML_COSB"][1].value,interval["D"][0].value,interval["D"][1].value)
    ax.text(xpos, ypos, words, size=.75*my_fontsize, fontname="monospace")
    ax.set_axisbelow(True)
    plt.tight_layout()
    # save the plot 
    outname = "../galactic_coordinate_plots/"+stream_name+"_"+cluster+"_"+cbar_name+".png"
    fig.savefig(outname, dpi=200)

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


# USE CHANGE THESE PARAMETERS
stream_name     = "gaia-11"
criterian_name  = "time_limits"
quantity_name   = "PML_COSB_short_cbar"
# PML_COSB_short_cbar, PML_COSB_long_cbar, PMB_short_cbar, PMB_long_cbar, time, distance

# set path information to the orbits of each cluster
back_path       = "../the_orbits/backward"
forward_path    = "../the_orbits/forward"
output_path     = "../outputs/" 
# load the interval
interval_fname  = "../outputs/"+stream_name+".json"
fp                  = open(interval_fname, "r+")
interval            = json.load(fp)
# generic plot properties
dot_size = .5

if "matches" in interval.keys():
    print("DOING {}".format(stream_name))
    # cluster's we're going to use if any are found within the range
    clusters = list(interval['matches'].keys())    
    # begin plot
    my_fontsize = 20
    fig,ax = plt.subplots(1,1)
    fig.set_figheight(5)
    fig.set_figwidth(11)
    # create the shape showing the latitude and longitude range
    lat_width = interval['LAT'][1] - interval['LAT'][0]
    lon_width = interval['LONG'][1] - interval['LONG'][0]
    me_patch=patches.Rectangle((interval['LONG'][0], interval['LAT'][0]), lon_width, lat_width, alpha=.4)
    # make plot nice
    ax.add_patch(me_patch)
    ax.set_xlim(-180,180)
    ax.set_ylim(-90,90)
    ax.set_xticks(np.arange(-180,210,30))
    ax.set_yticks(np.arange(-90,120,30))
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_aspect('equal')
    ax.set_title(stream_name, size=my_fontsize)
    ax.set_xlabel("LONG", size=my_fontsize)
    ax.set_ylabel("LAT",size=my_fontsize)
    ax.tick_params(labelsize=.75*my_fontsize)
    # tell the users of the criterai used in this region
    xpos = 180 - .05*(360)
    ypos = -90 + .02*(180)
    words = "PMB: [{:.1f} {:.1f}] mas/yr\nPML: [{:.1f} {:.1f}] mas/yr\nD:   [{:.1f} {:.1f}] kpc".format(interval['PMB'][0],interval["PMB"][1], interval["PML_COSB"][0], interval["PML_COSB"][1],interval["D"][0],interval["D"][1])
    ax.text(xpos, ypos, words, size=.75*my_fontsize, fontname="monospace")    
    if (type(interval["matches"]) is dict):
        # cluster's we're going to use if any are found within the range
        clusters = list(interval['matches'].keys())   
        for i in range(len(clusters)):
            # get the time and positions of each orbit
            t,x,y,z,vx,vy,vz = my_gf.concatenate_back_and_forward_orbit(clusters[i],back_path,forward_path)
            # put object into skycoord object
            galacto_centric = coords.Galactocentric(x=x,y=y,z=z,v_x=vx,v_y=vy, v_z=vz)
            # transform to galactic coordinates
            galactic_coordinates = galacto_centric.transform_to(Galactic)
            my_lon_coords = galactic_coordinates.l.wrap_at(180*u.degree)
            # set some generic properties
            my_cmap = "jet_r"
            # get the positions that passed the criterian test
            net_criterian = np.array([time in interval['matches'][clusters[i]] for time in t])
            # get the criterian to connect the cluster today to where its orbit was good 
            if criterian_name =="time_limits":
                criterian  = (interval["matches"][clusters[i]][0]  < t) & (t   < interval["matches"][clusters[i]][-1])
            elif criterian_name =="smart":
                connect_to_the_past = np.zeros(t.shape, dtype=bool)
                connect_to_the_future=np.zeros(t.shape, dtype=bool)
                the_time_of_my_life = np.array(interval['matches']['NGC5272'])
                # connect to the minimum negative value
                negative_time = the_time_of_my_life[the_time_of_my_life < 0]
                connect_to_the_past = (np.max(negative_time) < t) & (t < 0)
                # connect to the minimum positive value
                positive_time = the_time_of_my_life[the_time_of_my_life > 0]
                connect_to_the_future= (0 < t) & (t < np.min(positive_time))
                criterian = connect_to_the_future | connect_to_the_past
                criterian =  connect_to_the_future
            if quantity_name=="PML_COSB_long_cbar":
                cbar_name = "PML_COSB (mas/yr)"
                vmin = -20
                vmax =  20
                my_quantity         = galactic_coordinates.pm_l_cosb[criterian]
                good_quantity       = galactic_coordinates.pm_l_cosb[net_criterian]
            elif quantity_name== "PML_COSB_short_cbar":
                cbar_name = "PML_COSB (mas/yr)"
                vmin = -10
                vmax =  10
                my_quantity  = galactic_coordinates.pm_l_cosb[criterian]
                good_quantity= galactic_coordinates.pm_l_cosb[net_criterian]
            elif quantity_name == "PMB_short_cbar":
                cbar_name = "PMB (mas/yr)"
                vmin = -10
                vmax =  10
                my_quantity         = galactic_coordinates.pm_b[criterian]
                good_quantity         = galactic_coordinates.pm_b[net_criterian]
            elif quantity_name == "PMB_long_cbar":
                cbar_name = "PMB (mas/yr)"
                vmin = -20
                vmax =  20
                my_quantity         = galactic_coordinates.pm_b[criterian]
                good_quantity       = galactic_coordinates.pm_b[net_criterian]
            elif quantity_name =="time":
                cbar_name = "time (Gyrs)"
                my_quantity = t[criterian]
                good_quantity = t[net_criterian]
                vmin = -1
                vmax =  1
            elif quantity_name =="distance":
                my_quantity = galactic_coordinates.distance[criterian]
                good_quantity = galactic_coordinates.distance[net_criterian]
                vmin = 1.5
                vmax = 20
                my_cmap     = "viridis"
                cbar_name   = "Heliocentric distance (kpc)"
            # get coordinates
            my_x_coordinate     = my_lon_coords[criterian]
            my_y_coordinate     = galactic_coordinates.b[criterian]
            good_x_coordinates  = my_lon_coords[net_criterian]
            good_y_coordinates  = galactic_coordinates.b[net_criterian]            
            # plot the orbit of our globular cluster
            im = ax.scatter(my_x_coordinate, my_y_coordinate, c=my_quantity, s=dot_size, zorder=2,cmap=my_cmap, vmin=vmin, vmax=vmax)
            im_good = ax.scatter(good_x_coordinates,good_y_coordinates,c=good_quantity,s=20*dot_size,zorder=3,cmap=my_cmap, vmin=vmin, vmax=vmax)
            if i == 0:
                cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
                cbar.ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
                cbar.ax.get_yaxis().labelpad = 20
            # plot the current day position of our cluster
            index_cluster_today = np.argmin(abs(t))
            plt.scatter(my_lon_coords[index_cluster_today],galactic_coordinates.b[index_cluster_today], marker='o', c='k', zorder=4, alpha=.5)
            # add the name of the cluster
            plt.text(my_lon_coords[index_cluster_today].value+1,galactic_coordinates.b[index_cluster_today].value,clusters[i],size=.75*my_fontsize,fontname="monospace") 
    # save the figure
    ax.set_axisbelow(True)
    plt.tight_layout()    
    outname = "../stream_matches/"+stream_name+"_"+quantity_name+".png"
    fig.savefig(outname, dpi=200)
else:
    print("No matches reported for {}.\nRun intersection_json.py for this cluster".format(stream_name))
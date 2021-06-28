import numpy as np 
import matplotlib.pyplot as plt 
import h5py
import matplotlib.patches as patches
import json
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u 
from matplotlib.colors import LogNorm
import sys
from PIL import Image, ImageOps, ImageDraw, ImageFont
import os 
sys.path.append('../../functions')
import global_functions as my_gf 

# function the re-extracts the intersection of each list 
def extract_intersection(t, my_lon_coords, galactic_coordinates, interval):
    """ find which points we will color bold """
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
    return net_criterian

def do_plot(stream_name, cluster, interval, orbital_criterian, constrain_criterian,\
    my_lon_coords, galactic_coordinates, GDR2=True, dot_size=0.5):
    """ set the generic plot properties """
    # Change to match ibata et al depending on the stream 
    if GDR2==True:
        pm_vmin = -20
        pm_vmax = 20
    else:
        pm_vmin = -10
        pm_vmax =  10
    my_fontsize = 20
    fig,ax = plt.subplots(3,1)
    fig.set_figheight(15)
    fig.set_figwidth(11)
    # create the shape showing the latitude and longitude range
    lat_width = interval['LAT'][1] - interval['LAT'][0]
    lon_width = interval['LONG'][1] - interval['LONG'][0]
    me_patch= []
    # make plot nice
    for i in range(3):
        me_patch.append(patches.Rectangle((
            interval['LONG'][0].value, interval['LAT'][0].value),
            lon_width.value, lat_width.value, alpha=.4))
        ax[i].add_patch(me_patch[i])
        ax[i].set_xlim(-180,180)
        ax[i].set_ylim(-90,90)
        ax[i].set_yticks(np.arange(-90,120,30))
        ax[i].grid(True)
        ax[i].invert_xaxis()
        ax[i].set_aspect('equal')
        ax[i].set_ylabel("LAT",size=my_fontsize)
        ax[i].set_xticks(np.arange(-180,210,30))
        ax[i].tick_params(labelsize=.75*my_fontsize)
    ax[0].axes.xaxis.set_ticklabels([])
    ax[1].axes.xaxis.set_ticklabels([])
    ax[0].set_title("{:s}-{:s}".format(stream_name,cluster), size=my_fontsize)
    ax[2].set_xlabel("LONG", size=my_fontsize)
    # where to put the words of constraints used
    xpos = 180 - .05*(360)
    ypos = -90 + .02*(180)
    words = "PMB: [{:.1f} {:.1f}] mas/yr\nPML: [{:.1f} {:.1f}] mas/yr\nD:   [{:.1f} {:.1f}] kpc".format(\
        interval['PMB'][0].value,interval["PMB"][1].value, interval["PML_COSB"][0].value, \
        interval["PML_COSB"][1].value,interval["D"][0].value / 1000,interval["D"][1].value / 1000)
    ax[2].text(xpos, ypos, words, size=.75*my_fontsize, fontname="monospace")
    # plot the data each
    orbital_x     = my_lon_coords[orbital_criterian]
    orbital_y     = galactic_coordinates.b[orbital_criterian]
    constraint_x  = my_lon_coords[constrain_criterian]
    constraint_y  = galactic_coordinates.b[constrain_criterian]  
    im = []
    cbar = [] 
    for i in range(3): 
        if i == 2:
            main_orbit  = galactic_coordinates.pm_b[orbital_criterian]
            constraint  = galactic_coordinates.pm_b[constrain_criterian]
            my_cmap = "jet_r"
            vmin = pm_vmin
            vmax = pm_vmax
            cbar_name = "PMB (mas/yr)"
        elif i == 1:
            main_orbit  = galactic_coordinates.pm_l_cosb[orbital_criterian]
            constraint  = galactic_coordinates.pm_l_cosb[constrain_criterian]
            my_cmap = "jet_r"
            vmin = pm_vmin
            vmax = pm_vmax
            cbar_name = "PML cos(b) (mas/yr)"
        elif i ==0:
            main_orbit  = galactic_coordinates.distance[orbital_criterian]
            constraint  = galactic_coordinates.distance[constrain_criterian]
            my_cmap = "viridis"
            vmin = 1.5
            vmax = 20
            cbar_name = "Heliocentric Distance (kpc)"
        # do the plot of each 
        im.append(ax[i].scatter(orbital_x, orbital_y, c=main_orbit, s=dot_size, \
            zorder=2,cmap=my_cmap, vmin=vmin, vmax=vmax))
        if np.sum(constrain_criterian) > 0:
            ax[i].scatter(constraint_x,constraint_y,c=constraint,s=20*dot_size,\
                zorder=3,cmap=my_cmap, vmin=vmin, vmax=vmax)            
        # add the colorbar for each 
        cbar.append(fig.colorbar(im[i], ax=ax[i],fraction=0.046, pad=0.04))
        cbar[i].ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
        cbar[i].ax.get_yaxis().labelpad = 20
    return fig

def save_figure(stream_name,cluster, out_directory, figure, GDR2=True):
    """ save the figure """
    # check if directory exists
    if not os.path.isdir(out_directory+stream_name):
        os.mkdir(out_directory+stream_name+"/")
    figure.tight_layout()
    figure.savefig(out_directory+stream_name+"/"+cluster+".png", dpi=300)

def main():
    # function that opens 
    galactic_potential = "PII"
    back_path       = "../the_orbits/backward"
    forward_path    = "../the_orbits/forward"
    stream_name     = "fimbulthul" 
    interval_fname  = "../outputs/"+stream_name+".json"
    out_directory   = "../output_stream_to_cluster/"
    fp                  = open(interval_fname, "r+")
    interval            = json.load(fp)
    # cheeck to see if there are any matches
    if "matches" in interval.keys():
        # now see if there are any matches 
        if (type(interval["matches"]) is dict):
            matched_clusters = list(interval['matches'].keys())
            # now only get the clusters where we have the data files 
            common_clusters, _nothing = my_gf.match_cluster_list_to_data_directory(\
                matched_clusters,back_path,stream_or_backwardorbit="backwardorbits")
            # put our inteval object into the proper units 
            interval['LONG']    *=u.degree
            interval['LAT']     *= u.degree
            interval['PMB']     *=u.mas / u.year
            interval['PML_COSB']*=u.mas / u.year
            interval['D']       *=1000*u.pc
            for cluster in common_clusters:
                print(cluster)
                # do a plot for each cluster and save it 
                t,x,y,z,vx,vy,vz = my_gf.concatenate_back_and_forward_orbit(\
                    cluster, back_path, forward_path)
                galacto_centric = coord.Galactocentric(x=x,y=y,z=z,v_x=vx,v_y=vy, v_z=vz)
                # transform to galactic coordinates
                galactic_coordinates    = galacto_centric.transform_to(coord.Galactic)
                my_lon_coords           = galactic_coordinates.l.wrap_at(180*u.degree)
                constraint_criteran     = extract_intersection(\
                    t, my_lon_coords, galactic_coordinates, interval)
                all_points_criteran     = (interval["time"][0]  < t) &\
                    (t   < interval["time"][1])
                fig = do_plot(stream_name, cluster, interval, all_points_criteran, constraint_criteran,\
                    my_lon_coords, galactic_coordinates, GDR2=True, dot_size=0.5)
                print("Saving {:s}-{:s}".format(stream_name, cluster))
                save_figure(stream_name,cluster,out_directory,fig)
                
        else:
            print("No matches reported for {:s}".format(stream_name))
    else:
        print("Run Matches on {:s}".format(stream_name))    

if __name__ == "__main__":
    main()
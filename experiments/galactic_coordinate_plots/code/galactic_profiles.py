import numpy as np 
import matplotlib.pyplot as plt 
import h5py
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u 
import os 
import re
import sys
sys.path.append('../../functions')
import global_functions as my_gf 

# Get all of the streams
def get_stream(simulation_path, cluster_name, n_years_ago):
    """ Get the position, time, and stuff of the stream """ 
    data_file = h5py.File(simulation_path+cluster_name+".h5", 'r')
    criterian = data_file['energy']['Tescape'][:]/10 > 5-n_years_ago 
    lon = data_file['galactic']['LONG'][criterian]
    lat = data_file['galactic']['LAT'][criterian]
    dist= data_file['equatorial']['D'][criterian]
    return lon, lat, dist

def get_COM_coordinates(common_names):
    ''' take the clusters found in the directory and get COMS from vasiliev '''
    vasiliev=fits.open("../../simulations/cluster_data/asu_with_baumgardt_names.fit")
    # only extract the ones we are going to use        
    names = vasiliev[1].data['Name'].replace(' ', '') 
    RA = []
    dec = []
    dist = []
    # coordinate between the two lists
    for gc_name in common_names:
        for i, j in enumerate(names):
            if j == gc_name:
                RA.append(vasiliev[1].data['RAJ2000'][i])
                dec.append(vasiliev[1].data['DEJ2000'][i])
                dist.append(vasiliev[1].data['Dist'][i])
    # store in a cluster object
    glob_center =coord.SkyCoord(
        ra      = RA   *   u.degree,
        dec     = dec  *   u.degree, 
        distance= dist *   u.kpc,
        frame   = 'icrs')   
    # put this into galactic coordinates 
    galactic = glob_center.transform_to('galactic')
    cluster_lon = galactic.l
    cluster_lon[cluster_lon > 180*u.degree] = cluster_lon[cluster_lon > 180*u.degree] - 360*u.degree
    cluster_lat = galactic.b 
    return cluster_lon.value, cluster_lat.value, galactic.distance.value


def heliocentric_slice_plot(galactic_lon, galactic_lat, helio_centric_distance,
    down_distance, up_distance, 
    cluster_lon, cluster_lat, cluster_distance, cluster_names, outbasename):
    """ make the nice plot the way we love it """
    # set the path information
    outpath = "../output_heliocentric_profile/"
    outname = outbasename+"_heliocentric_distance_{}_{}_kpc.png".format(int(np.floor(down_distance)),int(np.ceil(up_distance)))
    cbar_name = "Heliocentric Distance (kpc)"
    # limits for the axis
    xmin = -180
    xmax =  180
    ymin = -90
    ymax =  90
    vmax = 30
    vmin = 0
    # other plot properties
    my_fontsize=20
    dot_size = .1
    fig, ax = plt.subplots(1,1)
    fig.set_figheight(5)
    fig.set_figwidth(11)
    # plot the data
    im = ax.scatter(galactic_lon, galactic_lat, c=helio_centric_distance, s=dot_size, vmin=vmin, vmax=vmax)
    # add the colorbar
    cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
    cbar.ax.get_yaxis().labelpad = 20
    # make the axes nice
    ax.grid(True)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    # make consistent with ibata 2021
    ax.invert_xaxis()
    ax.set_xticks(np.arange(-180,210,30))
    ax.set_yticks(np.arange(-90,120,30))
    # add the plot information 
    ax.set_xlabel("Longitude",   size=my_fontsize)
    ax.set_ylabel("Latitude",   size=my_fontsize)
    for i in range(len(cluster_names)):
        if (down_distance<cluster_distance[i]) & (cluster_distance[i] < up_distance):
            if cluster_lon[i] > 180:
                my_lon = cluster_lon[i]  - 360
            else:
                my_lon = cluster_lon[i]
            ax.text(my_lon, cluster_lat[i], cluster_names[i])
    ax.set_title("[{}, {}] kpc".format(down_distance, up_distance), size=my_fontsize)
    ax.tick_params(labelsize=.75*my_fontsize)
    ax.set_axisbelow(True)
    plt.tight_layout()
    fig.savefig(outpath+outname, dpi=150)
    plt.close()

def main():
    cluster_path = "../../functions/"
    clusters_fname = "clean_streams"
    simulation_path = "../../simulations/outputDATA/PII/streams/"
    common_clusters, n_years_ago=my_gf.match_cluster_list_to_data_directory(cluster_path+clusters_fname+".txt", \
        simulation_path, stream_or_backwardorbit="streams")
    # initialize the slicing spacing
    down_lim_list   = np.arange(0,25,1)
    up_lim_list     = np.arange(5,30,1)
    down_lim_list   = np.append(down_lim_list, 25)
    up_lim_list     = np.append(up_lim_list,   300)
    # initialize the lists for the stars
    lon = []
    lat = []
    dist = []
    # put the positions of the stars together
    for i in range(len(common_clusters)):
        my_lon, my_lat, my_dist = get_stream(simulation_path, common_clusters[i], n_years_ago[i])
        lon.append(my_lon)
        lat.append(my_lat)
        dist.append(my_dist)
    # turn list of lists into one numpy array
    longitudes  = np.hstack(lon)
    latitudes   = np.hstack(lat)
    distances   = np.hstack(dist)
    # get the name of each globular cluster and its distance from the center
    cluster_lon, cluster_lat, cluster_distance=get_COM_coordinates(common_clusters)
    # do each plot 
    for i in range(len(down_lim_list)):
        critera = (down_lim_list[i] < distances) &  (distances < up_lim_list[i])
        galactic_lon = longitudes[critera]
        galactic_lat = latitudes[critera]
        helio_centric_distance = distances[critera]
        heliocentric_slice_plot(galactic_lon, galactic_lat, helio_centric_distance, \
            down_lim_list[i], up_lim_list[i],\
                cluster_lon, cluster_lat, cluster_distance, common_clusters, clusters_fname)
    


if __name__=="__main__":
    main()

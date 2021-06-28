import numpy as np 
import h5py
import inspect
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u 
import os 
import re
import sys

def match_cluster_list_to_data_directory(cluster_list_or_path,simulation_path,stream_or_backwardorbit="streams"):
    """ 
    PURPOSE:
        Given a two column text file of cluster names and number of billion years, 
        match with the directory of .h5 and return the matches
    ARGUMENTS:
        cluster_list_or_path:
            path/file name of the text file of the clusters
            OR it is a python list of the clusters 
        simulation_path:
            Path to the folder containing data files of the clusters, i.e. backward orbit  
        stream_or_backwardorbit:
            Optional parameter that asks for the type of files. 
    """
    # Load in the cluster names 
    if type(cluster_list_or_path)==list:
        my_cluster_names=np.array(cluster_list_or_path)
        n_years_ago_list = 5*np.ones((len(cluster_list_or_path)))
    else:
        big_list            = np.loadtxt(cluster_list_or_path, dtype='str')
        my_cluster_names    = np.array([x[0] for x in big_list])            # names of the clusters
        n_years_ago_list    = np.array([int(x[1]) for x in big_list])       # number of years ago
    # check what our file type is
    if stream_or_backwardorbit=="streams":
        extension = ".*h5"
    elif stream_or_backwardorbit=="backwardorbits":
        extension = ".*dat"
    else:
        print("Warning, optional argument for 'stream_or_backwardorbit' is invalid\n {} is not 'streams' or 'backwardorbits'".format(stream_or_backwardorbit))
    # see if the path works and then get all of the clusters in the directory 
    if os.path.isdir(simulation_path):
        mylist  = os.listdir(simulation_path)
        r       = re.compile(extension)
        fnames  = list(filter(r.match, mylist))
        if stream_or_backwardorbit=="backwardorbits":
            fnames  = np.array([x.split('.')[0][5:] for x in fnames]) # skip the "orbit" in the file name
        else:
            fnames  = np.array([x.split('.')[0] for x in fnames])
    else:
        print("{} was not found ".format(simulation_path))
    # find the intersection set the of the files in the directory
    common_elements = np.in1d(fnames, my_cluster_names)
    new_list=fnames[common_elements]
    common_n_years_ago      =   n_years_ago_list[np.in1d(my_cluster_names, new_list)]
    common_cluster_names    =   my_cluster_names[np.in1d(my_cluster_names, new_list)]
    # return the common elements and the number of years ago in the correct order
    return common_cluster_names, common_n_years_ago


def concatenate_back_and_forward_orbit(cluster, back_path, forward_path):
    ''' 
    Purpose:
        Concatenate the backward and forward orbit together
    Arguments:
        cluster:        STRING. the name of the cluster
        back_path:      STRING. the path to the folder with orbits back in time
        forward_path:   STRING. the path to the folder with the orbits foward in time
    Returns:
        t, x, y, z, vx, vy, vz 
        t is in bilions of years
        x,y,z are in kpc
        vx,vy,vz are in km/s
    # open forward and backward orbits
    '''
    print("forward_path", forward_path)
    print(cluster)
    print(forward_path+"/"+"orbit"+cluster+".dat")
    tf,xf,yf,zf,vxf,vyf,vzf = np.loadtxt(forward_path+"/"+"orbit"+cluster+".dat", unpack=True,skiprows=1)
    tb,xb,yb,zb,vxb,vyb,vzb = np.loadtxt(back_path+"/"+"orbit"+cluster+".dat", unpack=True,skiprows=1)
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
    return t,x,y,z,vx,vy,vz 

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

# Get all of the streams
def get_stream(simulation_path, cluster_name, n_years_ago=5):
    """ Report the positions and velocities of all the escapes stars 
    Arguments:
    simulation_path : the path to the folder containig the streams
    cluster_name    : the name of the cluster, without the .h5 extension
    n_years_age = 5 : time filter for stream in billions of years ago. Defaults to 5, 
                        which allows all stars to enter
    RETURNS
    x,y,z       = kpc
    vx,vy,vz    = 10 km/s 
    """ 
    data_file = h5py.File(simulation_path+cluster_name+".h5", 'r')
    criterian = data_file['energy']['Tescape'][:]/10 > 5-n_years_ago 
    x = data_file['galactocentric']['X'][criterian]
    y = data_file['galactocentric']['Y'][criterian]
    z = data_file['galactocentric']['Z'][criterian]
    vx = data_file['galactocentric']['VX'][criterian]
    vy = data_file['galactocentric']['VY'][criterian]
    vz = data_file['galactocentric']['VZ'][criterian]
    return x,y,z,vx,vy,vz

def get_backward_orbit(path):
    """ provide the path to the cluster's .dat file of its orbit """
    # the orbit of the body in time
    t,x,y,z,vx,vy,vz = np.loadtxt(path+".dat", unpack=True)
    t = np.flip(t/10)  # Now in billions of years ago
    
    x = np.flip(x)  
    y = np.flip(y)
    z = np.flip(z)
    vx = np.flip(vx)
    vy = np.flip(vy)
    vz = np.flip(vz)
    return t, x, y, z, vx, vy, vz 


def get_COM_coordinates_RA(common_names):
    ''' take the clusters found in the directory and get COMS from vasiliev '''
    vasiliev=fits.open("../../simulations/cluster_data/asu_with_baumgardt_names.fit")
    # only extract the ones we are going to use        
    names = vasiliev[1].data['Name'].replace(' ', '') 
    RA = []
    dec = []
    dist = []
    found_names = []
    # coordinate between the two lists
    for gc_name in common_names:
        for i, j in enumerate(names):
            if j == gc_name:
                RA.append(vasiliev[1].data['RAJ2000'][i])
                dec.append(vasiliev[1].data['DEJ2000'][i])
                dist.append(vasiliev[1].data['Dist'][i])
                found_names.append(gc_name)
    # store in a cluster object
    glob_center =coord.SkyCoord(
        ra      = RA   *   u.degree,
        dec     = dec  *   u.degree, 
        distance= dist *   u.kpc,
        frame   = 'icrs')  
    return glob_center, found_names
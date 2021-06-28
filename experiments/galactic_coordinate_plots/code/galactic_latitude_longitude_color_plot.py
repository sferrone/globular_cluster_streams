import numpy as np 
import matplotlib.pyplot as plt 
import h5py
# import astropy.coordinates as coord
import os 
import re

# USER can changes this information 
cluster_list        = "mini_list_of_clusters.txt.txt"    # file name of strings of desired clusters
quantity_name       = 'PMB'                     # data quantity to plot 
galactic_potential  = "PII"                     # galactic potential model 
globular_potential  = "Plummer"                 # globular cluster potential
dot_size = 0.2                                  # size of dot for scatter plot 
quantity_coordinate_system = 'galactic'         # units from h5 data table
coordinate_system          = 'galactic'         # units for the axes
# how many years ago are we going to go?
recent=1
n_years_ago=1
'''
equatorial
['D', 'DEC', 'PMDEC', 'PMRA_COSDEC', 'RA']
galactic
['LAT', 'LONG', 'PMB', 'PML_COSB']
galactocentric
['VX', 'VY', 'VZ', 'X', 'Y', 'Z']
^ Also Lz
'''
# initialize quantities 
simulation_path = "../../simulations/outputDATA/"+galactic_potential+"/"+globular_potential+"/streams/"
outpath = "../outputs/"
title_name = r"{0}-{1}-$\bf{2}$".format(galactic_potential,globular_potential, cluster_name)
position_sun_galactic_coordaintes = np.array([-8.34,0,0.027])
# initialie master list for all
quantity    = []
xcoords     = []
ycoords     = []
cluster_xcoords =   []
cluster_ycoords =   []
# Get the list of streams for each cluter in our directory
mylist  = os.listdir(simulation_path)
r       = re.compile(".*h5")
fnames  = list(filter(r.match, mylist))
fnames  = [x.split('.')[0] for x in fnames]
# get clusters from text file
fp = open(cluster_file_name)
my_cluster_names=fp.readlines()
my_cluster_names = [x.split('\n')[0] for x in my_cluster_names]
# find the intersection set of desired clusters and those in the directory
common_elements = list(set(my_cluster_names) & set(fnames))
n_files = int(len(common_elements))
# open each h5 file and take the wanted data
for cluster_name in common_elements:
    data_file = h5py.File(simulation_path+cluster_name, 'r')
    if recent==1:
        criterian = data_file['energy']['Tescape'][:]/10 > 5-n_years_ago
        print("{} stars within the last {} billion years ".format(np.sum(criterian), n_years_ago))
    else:
        criterian = np.ones(data_file['energy']['Tescape'][:].shape, dtype=bool)
    # extract the RA and DEC of the cluster
    vasiliev=fits.open("../../simulations/cluster_data/asu.fit")
    alpha_GC=vasiliev[1].data['RAJ2000']
    delta_GC=vasiliev[1].data['DEJ2000']
    D       =vasiliev[1].data['Dist']
    # get index of cluster in vasilief
    for i,j in enumerate(vasiliev[1].data['Name'].replace(' ','')):
        if j == cluster_name:
            my_index = i
    # create the coordinate object        
    glob_center =coord.SkyCoord(
        ra      = alpha_GC[my_index] *   u.degree,
        dec     = delta_GC[my_index] *   u.degree, 
        distance= D[my_index]        *   u.kpc,
        frame   = 'icrs')       
    # decide which axes to plot based on the coordinate system
    if coordinate_system=="galactic":
        xlabel = "LONG"
        ylabel = "LAT"
        xmin = -180
        xmax =  180
        ymin = -90
        ymax =  90    
        my_x_coordinate = np.array(data_file[coordinate_system][xlabel][criterian])
        my_y_coordinate = np.array(data_file[coordinate_system][ylabel][criterian])
        galactic = glob_center.transform_to('galactic')
        cluster_X = galactic.l
        cluster_Y = galactic.b
    elif coordinate_system=="equatorial":    
        xlabel = "RA"
        ylabel = "DEC" 
        xmin =  0
        xmax =  360
        ymin = -90
        ymax =  90
        # store positions of all stars
        my_x_coordinate = np.array(data_file[coordinate_system][xlabel][criterian])
        my_y_coordinate = np.array(data_file[coordinate_system][ylabel][criterian])
        ### Get the cluster's center position
        cluster_X = glob_center.ra
        cluster_Y = glob_center.dec
    else:
        xlabel = "X (KPC)"
        ylabel = "Y (KPC)" 
        my_x_coordinate = np.array(data_file[coordinate_system]["X"][criterian])
        my_y_coordinate = np.array(data_file[coordinate_system]["Y"][criterian])    
        # pick limits by rounding up to the nearest 10   
        extremes = [np.max(abs(my_x_coordinate)), np.max(abs(my_y_coordinate))]
        my_extreme = round(np.max(extremes)/10)*10
        xmin = -my_extreme
        xmax =  my_extreme
        ymin = -my_extreme
        ymax =  my_extreme
        # get the position in galacto-centric distance
        gc_frame = ref_frame()
        galactocentric = glob_center.transform_to(gc_frame) 
        cluster_X = galactocentric.x
        cluster_Y = galactocentric.y
    # Get colorbar based on the unit selected 
    my_cmap = "jet_r"
    strict_color_range = 0
    cbar_name = quantity_name
    # calculate or extract which quantity to use as color scale
    if quantity_name == "Lz":
        X = np.array(data_file['galactocentric']['X'][criterian])
        Y = np.array(data_file['galactocentric']['Y'][criterian])
        VX = np.array(data_file['galactocentric']['VX'][criterian])
        VY = np.array(data_file['galactocentric']['VY'][criterian])
        angular_momentum = X*VY - Y*VX
        my_quantity = np.array(angular_momentum)
        cbar_name = "Z-angular momentum (kpc*kpc/year)"
    elif quantity_name == "D":
        cbar_name = "Heliocentric distance (kpc)"
        vmin = 1.5
        vmax = 20
        strict_color_range = 1
        my_cmap = "viridis"
        my_quantity = np.array(data_file[quantity_coordinate_system][quantity_name][criterian])
    elif quantity_name == "PMB":
        vmin = -20
        vmax =  20
        strict_color_range = 1
    elif quantity_name == "PML_COSB":    
        vmin = -20
        vmax =  20
        strict_color_range = 1    
    elif quantity_name=="t_escape":
        cbar_name = "Time-Escape (Gyrs ago)"
        criterian               =   data_file['energy']['Tescape'][:]/10 > (5-n_years_ago)
        time_escape_fugitives   =   data_file['energy']['Tescape'][criterian]/10 # put in billions of years
        my_x_coordinate = np.array(data_file[coordinate_system][xlabel][criterian])
        my_y_coordinate = np.array(data_file[coordinate_system][ylabel][criterian])
        my_x_coordinate_grey    = np.array(data_file[coordinate_system][xlabel][~criterian])
        my_y_coordinate_grey    = np.array(data_file[coordinate_system][ylabel][~criterian])
        my_quantity             = time_escape_fugitives - 5    # put this in billion of years ago       
        vmax=0
        vmin=-n_years_ago
        strict_color_range=1
    else:
        my_quantity = np.array(data_file[quantity_coordinate_system][quantity_name][criterian])
    quantity.append(my_quantity)
    xcoords.append(my_x_coordinate)
    ycoords.append(my_y_coordinate)
    cluster_xcoords.append(cluster_X)
    cluster_ycoords.append(cluster_Y)

# begin plot 
fig, ax = plt.subplots(1,1)
my_fontsize=20
fig.set_figheight(5)
fig.set_figwidth(11)
if coordinate_system=="galactocentric":
    fig.set_figheight(5)
    fig.set_figwidth(5)
else:
    fig.set_figheight(5)
    fig.set_figwidth(11)
ax.grid(True)
# plot the data 
ax.plot(cluster_X, cluster_Y, 'kx', zorder=3, markersize=8)
if strict_color_range==0:
    im = ax.scatter(my_x_coordinate, my_y_coordinate, c=my_quantity, s=dot_size, zorder=2,cmap=my_cmap)
elif strict_color_range==1:
    im = ax.scatter(my_x_coordinate, my_y_coordinate, c=my_quantity, s=dot_size, zorder=2,cmap=my_cmap, vmin=vmin, vmax=vmax)
# plot the grey points that have not escaped
if quantity_name=="t_escape":
    im_2 = ax.scatter(my_x_coordinate_grey, my_y_coordinate_grey, s=dot_size, zorder=1, color='k')
# make the plot nice
cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
cbar.ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
cbar.ax.get_yaxis().labelpad = 20
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect('equal')
# add information 
ax.set_title(title_name, size=my_fontsize)
ax.set_xlabel(xlabel, size=my_fontsize)
ax.set_ylabel(ylabel,size=my_fontsize)
ax.tick_params(labelsize=.75*my_fontsize)
ax.set_axisbelow(True)
plt.tight_layout()
# make consistent with rodrigo 2019
if coordinate_system=="galactic":
    ax.invert_xaxis()
    ax.set_xticks(np.arange(-180,210,30))
    ax.set_yticks(np.arange(-90,120,30))
# save plot and write out our data
n_files = len(fnames)
outname = coordinate_system+"_coordinates_"galactic_potential+"_"+globular_potential+"_"+quantity_name+"_nClusters_"+str(n_files)+".png"
fig.savefig(outpath+outname, dpi=150)

# close the data file
data_file.close()
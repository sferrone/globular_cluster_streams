import numpy as np 
import matplotlib.pyplot as plt 
import h5py
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u 
import os 
import re
import sys
# USER can changes this information 
cluster_file_name   = "mini_list_of_clusters.txt"   # file name of strings of desired clusters
quantity_name       = 'D'                           # data quantity to plot 
galactic_potential  = "PII"                         # galactic potential model 
globular_potential  = "Plummer"                     # globular cluster potential
dot_size = 0.2                                      # size of dot for scatter plot 
quantity_coordinate_system = 'equatorial'           # units from h5 data table
coordinate_system          = 'galactic'             # units for the axes
# how many years ago are we going to go?
recent=1
'''
equatorial
['D', 'DEC', 'PMDEC', 'PMRA_COSDEC', 'RA', 'RV']
galactic
['LAT', 'LONG', 'PMB', 'PML_COSB']
galactocentric
['VX', 'VY', 'VZ', 'X', 'Y', 'Z']
^ Also Lz
'''
# initialize quantities 
simulation_path = "../../simulations/outputDATA/"+galactic_potential+"/"+globular_potential+"/streams/"
outpath = "../outputs/"
title_name = r"{0}-{1}-$\bf{2}$".format(galactic_potential,globular_potential, "fun")
vasiliev=fits.open("../../simulations/cluster_data/asu.fit")
position_sun_galactic_coordaintes = np.array([-8.34,0,0.027])
# initialie master list for all
quantity    = []
xcoords     = []
ycoords     = []
cluster_xcoords =   []
cluster_ycoords =   []
geometric_mean_x=   []
geometric_mean_y=   []
# Get the list of streams for each cluter in our directory
mylist  = os.listdir(simulation_path)
r       = re.compile(".*h5")
fnames  = list(filter(r.match, mylist))
fnames  = np.array([x.split('.')[0] for x in fnames])
# get clusters from text file
big_list=np.loadtxt(cluster_file_name, dtype='str')
my_cluster_names = np.array([x[0] for x in big_list])
n_years_ago_list = np.array([int(x[1]) for x in big_list])
# find the intersection set of desired clusters and those in the directory
# common_elements = list(set(my_cluster_names) & set(fnames))
common_elements = np.in1d(fnames, my_cluster_names)
new_list=fnames[common_elements]
n_years_ago_list=n_years_ago_list[np.in1d(my_cluster_names, new_list)]
my_cluster_names=my_cluster_names[np.in1d(my_cluster_names, new_list)]
n_files = int(len(n_years_ago_list))
print(my_cluster_names)
# open each h5 file and take the wanted data
counter=0
for cluster_name in my_cluster_names:
    data_file = h5py.File(simulation_path+cluster_name+".h5", 'r')
    # filter data by number of years ago
    if recent==1:
        criterian = data_file['energy']['Tescape'][:]/10 > 5-n_years_ago_list[counter]
        print("{}: {} stars within the last {} billion years ".format(cluster_name, np.sum(criterian), n_years_ago_list[counter]))
    else:
        criterian = np.ones(data_file['energy']['Tescape'][:].shape, dtype=bool)
    # extract the RA and DEC of the cluster
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
        cluster_Y = galactic.b
        if galactic.l < 180*u.degree:
            cluster_X = galactic.l
        else:
            cluster_X =  galactic.l-360*u.degree
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
    counter+=1
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
        my_quantity = np.array(data_file[quantity_coordinate_system][quantity_name][criterian])
    elif quantity_name == "PML_COSB":    
        vmin = -20
        vmax =  20
        strict_color_range = 1   
        my_quantity = np.array(data_file[quantity_coordinate_system][quantity_name][criterian]) 
    elif quantity_name=="t_escape":
        cbar_name = "Time-Escape (Gyrs ago)"
        criterian               =   data_file['energy']['Tescape'][:]/10 > (5-n_years_ago_list[counter])
        time_escape_fugitives   =   data_file['energy']['Tescape'][criterian]/10 # put in billions of years
        my_x_coordinate = np.array(data_file[coordinate_system][xlabel][criterian])
        my_y_coordinate = np.array(data_file[coordinate_system][ylabel][criterian])
        my_x_coordinate_grey    = np.array(data_file[coordinate_system][xlabel][~criterian])
        my_y_coordinate_grey    = np.array(data_file[coordinate_system][ylabel][~criterian])
        my_quantity             = time_escape_fugitives - 5    # put this in billion of years ago       
        vmax=0
        vmin=-5
        strict_color_range=1
    else:
        my_quantity = np.array(data_file[quantity_coordinate_system][quantity_name][criterian])
    quantity.append(my_quantity)
    xcoords.append(my_x_coordinate)
    ycoords.append(my_y_coordinate)
    geometric_mean_x.append(np.mean(my_x_coordinate)*u.degree)
    geometric_mean_y.append(np.mean(my_y_coordinate)*u.degree)
    cluster_xcoords.append(cluster_X)
    cluster_ycoords.append(cluster_Y)
# put them into one diemsional list
quantity        = np.hstack(quantity)
xcoords         = np.hstack(xcoords)
ycoords         = np.hstack(ycoords)
cluster_xcoords = np.hstack(cluster_xcoords)
cluster_ycoords = np.hstack(cluster_ycoords)

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
ax.plot(cluster_xcoords, cluster_ycoords, 'kx', zorder=3, markersize=8)
if strict_color_range==0:
    im = ax.scatter(xcoords, ycoords, c=quantity, s=dot_size, zorder=2,cmap=my_cmap)
elif strict_color_range==1:
    im = ax.scatter(xcoords, ycoords, c=quantity, s=dot_size, zorder=2,cmap=my_cmap, vmin=vmin, vmax=vmax)
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
ax.set_title(title_name,size=my_fontsize)
ax.set_xlabel(xlabel,   size=my_fontsize)
ax.set_ylabel(ylabel,   size=my_fontsize)
ax.tick_params(labelsize=.75*my_fontsize)
ax.set_axisbelow(True)
plt.tight_layout()
# add the name of the cluster
for i in range(len(my_cluster_names)):
    GC_X = cluster_xcoords[i]
    GC_Y = cluster_ycoords[i]
    # do the x coordinate
    xpos = geometric_mean_x[i] #- 10*np.abs(geometric_mean_x[i]-GC_X)
    # do the y coordiante
    if geometric_mean_y[i] < GC_Y:
        ypos = geometric_mean_y[i] #- np.abs(geometric_mean_y[i]-GC_Y)
    else:
        ypos = geometric_mean_y[i] #+ np.abs(geometric_mean_y[i]-GC_Y)
    print(my_cluster_names[i],xpos, ypos)
    ax.text(xpos.value, ypos.value, my_cluster_names[i], size=12, zorder=4)
    
# make consistent with rodrigo 2019
if coordinate_system=="galactic":
    ax.invert_xaxis()
    ax.set_xticks(np.arange(-180,210,30))
    ax.set_yticks(np.arange(-90,120,30))
# save plot and write out our data
outname = coordinate_system+"_coordinates_"+galactic_potential+"_"+globular_potential+"_"+quantity_name+"_nClusters_"+str(n_files)+".png"
print(outname)
fig.savefig(outpath+outname, dpi=150)

# close the data file
data_file.close()
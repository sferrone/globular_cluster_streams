import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.stats import sigma_clip
import h5py

# import astropy.coordinates as coord
import os 
import re

# USER can changes this information 

quantity_name = 'Lz'                         # data quantity to plot 
galactic_potential = "PII"                  # galactic potential model 
globular_potential = "Plummer"              # globular cluster potential
dot_size = 0.2                              # size of dot for scatter plot 
quantity_coordinate_system = 'galactocentric'   # units from h5 data table
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
position_sun_galactic_coordaintes = np.array([-8.34,0,0.027])
my_quantity = []
my_x_coordinate = []
my_y_coordinate = []

# get the all of the files we want 
mylist=os.listdir(simulation_path)
r = re.compile(".*h5")
fnames=list(filter(r.match, mylist))

# open each h5 file and take the wanted data
for fname in fnames:
    data_file = h5py.File(simulation_path+fname, 'r')
    data_fields = list(data_file.keys())
    if quantity_name == "Lz":
        X = np.array(data_file['galactocentric']['X'])
        Y = np.array(data_file['galactocentric']['Y'])
        VX = np.array(data_file['galactocentric']['VX'])
        VY = np.array(data_file['galactocentric']['VY'])
        angular_momentum = X*VY - Y*VX
        my_quantity.append(angular_momentum)
    else:
        my_quantity.append(np.array(data_file[quantity_coordinate_system][quantity_name]))
        
    my_x_coordinate.append(np.array(data_file['galactocentric']['X']))
    my_y_coordinate.append(np.array(data_file['galactocentric']['Y']))
x_coordinate = np.hstack(my_x_coordinate)    
y_coordinate = np.hstack(my_y_coordinate)    
quantity = np.hstack(my_quantity)  

# Get colorbar based on the unit selected 
if quantity_name == "D":
    cbar_name = "Heliocentric distance (kpc)"
elif quantity_name == "Lz":
    cbar_name = "Z-angular momentum (kpc* 10 km/s)"
else:
    cbar_name = quantity_name

# begin plot 
fig, ax = plt.subplots(1,1)
my_fontsize=20
fig.set_figheight(7)
fig.set_figwidth(7)
ax.grid(True)
# plot the data 
im = ax.scatter(x_coordinate, y_coordinate, c=quantity, s=dot_size)
# make colorbar fit axes and also 
divider = make_axes_locatable(ax)
cax1 = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax1)
cbar.ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
cbar.ax.get_yaxis().labelpad = 20
ax.set_aspect('equal')
# do sig_clip for colorbar


ax.set_xlim(-75, 75)
ax.set_ylim(-75,75)
ax.set_title("Galacto-centric coordinates")
ax.set_xlabel("X (kpc)", size=my_fontsize)
ax.set_ylabel("Y (kpc)",size=my_fontsize)
ax.tick_params(labelsize=.75*my_fontsize)
ax.set_axisbelow(True)
plt.tight_layout()

# save plot and write out our data
n_files = len(fnames)
outname = "galactic_coordinates_colorplot_quantity_"+quantity_name+"_nClusters_"+str(n_files)+".png"
fig.savefig(outpath+outname, dpi=150)

# close the data file
data_file.close()
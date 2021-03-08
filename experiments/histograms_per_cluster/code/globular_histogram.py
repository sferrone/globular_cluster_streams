import numpy as np 
import matplotlib.pyplot as plt 
import h5py
import output
from astropy.io import fits

# USER INPUTS
coordinate_system = 'galactocentric'
quantity_name = 'Lz'
'''equatorial
['D', 'DEC', 'PMDEC', 'PMRA_COSDEC', 'RA']
galactic
['LAT', 'LONG', 'PMB', 'PML_COSB']
galactocentric
['VX', 'VY', 'VZ', 'X', 'Y', 'Z']
^ Also Lz
'''
# cluster and model info
galactic_potential = "PII"
globular_potential = "Plummer"
cluster_name = 'Pal12'

# define path variables and constants
simulation_path = "../../simulations/outputDATA/"+galactic_potential+"/"+globular_potential+"/streams/"
outpath = "../outputs"
position_sun_galactic_coordaintes = np.array([-8.34,0,0.027])

# OPEN THE DATA
data_file = h5py.File(simulation_path+cluster_name+".h5", 'r')
# choose which data to plot
if quantity_name == "Lz":
    X = np.array(data_file['galactocentric']['X'])
    Y = np.array(data_file['galactocentric']['Y'])
    VX = np.array(data_file['galactocentric']['VX'])
    VY = np.array(data_file['galactocentric']['VY'])
    my_data = X*VY - Y*VX
else:
    my_data = data_file[coordinate_system][quantity_name][:]

# make case based xlabel name
if quantity_name=="D":
    xlabel_name = "Helio-centric distance (kpc)"
elif quantity_name == "Lz":
    xlabel_name=r"$L_z$ angular momentum"
else:
    xlabel_name = quantity_name

fontsize = 15
fig,ax = plt.subplots(1,1)
my_n_bins = int(np.sqrt(my_data.shape[0]))
ax.hist(my_data, bins=my_n_bins);
ax.set_ylabel("Number of Stars", size=fontsize)
ax.set_xlabel(xlabel_name, size=fontsize)
ax.set_title("{}-{}-{}".format(galactic_potential, globular_potential, cluster_name) , size=fontsize)
ax.grid(True)

# store this file
out_path = "../outputs/"+galactic_potential+"/"+globular_potential+"/"+cluster_name+"/"
fname = "histogram_"+quantity_name+"_"+cluster_name+".png"
output.createdir(out_path)
fig.savefig(out_path+fname, dpi=300)

# close h5 file
data_file.close()
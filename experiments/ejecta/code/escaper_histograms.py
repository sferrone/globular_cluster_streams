import numpy as np 
import matplotlib.pyplot as plt 
import time
import h5py
import sys 
import os
import re

# get backward orbit
GCname              =   'Eridanus'
galactic_potential  =   'PII'
cluster_potential   =   'Plummer'
simulation_path     =   '../../simulations/outputDATA/'+galactic_potential+"/"+cluster_potential+"/"
path_backward       =   simulation_path+"backwardorbits/"+"orbit"+GCname+".dat"
path_streams        =   simulation_path+"streams/"+GCname+".h5"
out_name            =   '../outputs/escaper_histogram_'+galactic_potential+"_"+cluster_potential+"_"+GCname+".png"

# get the all of the files we want 
mylist=os.listdir(simulation_path+"streams/")
r = re.compile(".*h5")
fnames=list(filter(r.match, mylist))
# check to see if the item is within the list 
cluster_names = [x.split('.')[0] for x in fnames]

if GCname in cluster_names:
    print("Doing {0}".format(GCname))
    # the current stream positions 
    streams = h5py.File(path_streams, 'r')
    # the orbit of the body in time
    t,x,y,z,vx,vy,vz = np.loadtxt(path_backward, unpack=True)
    t = t/10  # Now in billions of years ago
    x = np.flip(x)  
    y = np.flip(y)
    r = np.sqrt(x**2 + y**2)
    R = np.sqrt(x**2 + y**2 + z**2)
    # # find out when the stars escape
    criterian               =   streams['energy']['Tescape'][:] > 0
    fugitives               =   np.where(criterian)[0]
    time_escape_fugitives   =   streams['energy']['Tescape'][fugitives]/10 # put in billions of years

    # get the position of the cluster when each star escapes
    radius_escapers = np.ones((time_escape_fugitives.shape))
    z_escapers      = np.ones((time_escape_fugitives.shape))
    R_escapers      = np.ones((time_escape_fugitives.shape))
    for i in range(time_escape_fugitives.shape[0]):
        mindex = np.argmin(abs(t - time_escape_fugitives[i]))
        radius_escapers[i]  = r[mindex]
        z_escapers[i]       = z[mindex]
        R_escapers[i]       = R[mindex]

    # star the plot
    fig, ax = plt.subplots(2,2, figsize=(10,5))
    n_bins = np.floor(np.sqrt(time_escape_fugitives.shape[0]));
    n_escap = int(np.sum(criterian))
    fontsize=15
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    fig.delaxes(ax[1,1])
    # plot the coordinates versus time
    ax[0,0].plot(t,R, color=(0,0,1))
    ax[0,0].plot(t,r, color=(1,0,0))
    ax[0,0].plot(t,z, color=(0,1,0))
    # make plot nice 
    ax[0,0].set_ylabel("Distance (kpc)", size=fontsize)
    ax[0,0].grid(True)

    # do the horizontal histogram of each coordinate
    hist1=ax[0,1].hist(R_escapers,bins=int(n_bins),     orientation="horizontal",lw=5, fc=(0, 0, 1, 0.5))
    hist2=ax[0,1].hist(radius_escapers,bins=int(n_bins),orientation="horizontal",lw=5, fc=(1, 0, 0, 0.5))
    hist3=ax[0,1].hist(z_escapers,bins=int(n_bins),     orientation="horizontal", lw=5, fc=(0, 1, 0, 0.5))
    # remove yticks and make graph nice
    ax[0,1].grid(True)
    ax[0,1].set_ylim(ax[0,0].get_ylim())
    ax[0,1].set_xticks(ax[0,1].get_xticks()[1::])
    ax[0,1].yaxis.set_ticklabels([])
    # add information
    ax[0,1].set_xlabel("# Escapers", size=fontsize)
    ax[0,1].legend(["Radius (Spherical)", "Radius (Cylindrical)", "Elevation (Z)"], prop={'size':fontsize},loc='upper center', bbox_to_anchor=(0.5, -0.4))
    # print the name of the cluster
    xmin,xmax = ax[0,1].get_xlim()
    ymin,ymax = ax[0,1].get_ylim()
    xpos = .15*(xmax-xmin) + xmin
    ypos = -.425*(ymax-ymin) + ymin
    ax[0,1].text(xpos,ypos, r"$\bf{{{0}}}$ total: {1}".format(GCname, n_escap), size=fontsize)
    
    # do the histogram in time 
    ax[1,0].hist(time_escape_fugitives, bins=int(n_bins));
    ax[1,0].set_ylabel("# Escapers", size=fontsize)
    ax[1,0].set_xlabel("Time (Gyrs)", size=fontsize)
    ax[1,0].grid(True)
    ax[1,0].set_xlim(ax[0,0].get_xlim())
    fig.tight_layout()
    fig.savefig(out_name, bbox_inches = "tight", dpi=150)

else:
    print("{0} not found. Skipping".format(GCname))
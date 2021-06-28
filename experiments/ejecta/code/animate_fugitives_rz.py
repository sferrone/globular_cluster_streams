import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import time 
import h5py
# import astropy.coordinates as coord
import os 
import re
import sys

# get backward orbit
GCname              = sys.argv[1]
galactic_potential  = 'PII'
cluster_potential   = 'Plummer'
path                = '../../simulations/outputDATA/'+galactic_potential+"/"+cluster_potential+"/"
path_backward       = path+"backwardorbits/"+"orbit"+GCname+".dat"
path_streams        = path+"streams/"+GCname+".h5"

outname = "../outputs/"+GCname+"_"+galactic_potential+"_"+cluster_potential+"_fugitives_rz_test.gif"
# the current stream positions 
streams = h5py.File(path_streams, 'r')
# the orbit of the body in time
t,x,y,z,vx,vy,vz = np.loadtxt(path_backward, unpack=True)

# create class object to store the information about the escapers
class Escapers:
    """ need to define class object in animation to keep track of new escapers """ 
    def __init__(self, time_escape):
        " initialize by saying no new stars have escaped "
        self.total_escapers = 0
        self.time_escape = time_escape # a numpy array escape times from .h5 stream files
        self.new_escapers = 0
    def get_new_escapers(self, input_time):
        " the number of new escapers only increases  "
        calculated_escapers = sum(self.time_escape < input_time)
        if (calculated_escapers > self.total_escapers):
            self.new_escapers = calculated_escapers - self.total_escapers
            self.total_escapers = calculated_escapers
        else:
            self.new_escapers = 0
    def reset_escapers(self):
        ''' maybe we will need to reset the object one day '''
        self.new_escapers = 0
        self.total_escapers = 0

# make time go from oldest to now
t = t/10  # Now in billions of years ago
# reverse each to make sure we're going in the correct order
x = np.flip(x)  
y = np.flip(y)
r = np.sqrt(x**2 + y**2)
z = np.flip(z)

# set parameters for which frames to use 
n_skip      = 1000
start_frame = 0
n_frames    = int(np.floor(len(x)))
my_frames   = range(start_frame, n_frames, n_skip)

# # find out when the stars escape
fugitives               = streams['energy']['Tescape'][:] > 0
time_escape_fugitives   = streams['energy']['Tescape'][fugitives]/10 # put in billions of years
# Create escapers object to handle when stars escape
my_escapers = Escapers(time_escape_fugitives)
# First set up the figure, the axis, and the plot element we want to animate
fig, ax     = plt.subplots(1,1)
fontsize    = 15
# set axes limits
z_max = 1.1*np.max(np.max(abs(z)))
z_min = -z_max
r_max = 1.1*np.max(r)
r_min = 0
ax.set_xlim(r_min, r_max)
ax.set_ylim(z_min ,z_max)
# get coordiantes for witing text
xpos = .675*(r_max-r_min) + r_min
ypos = .925*(z_max-z_min) + z_min
# make graph nice
ax.set_xlabel("Galactic R (kpc)")
ax.set_ylabel("Galactic Z (kpc)")
ax.set_aspect('equal')

# animation function.  This is called sequentially
def animate(i):
    # start plot
    ax.clear()
    # set plot properties
    ax.set_xlim(r_min, r_max)
    ax.set_ylim(z_min ,z_max)
    ax.set_xlabel("Galactic R (kpc)", size=fontsize)
    ax.set_ylabel("Galactic Z (kpc)", size=fontsize)    
    ax.set_title(GCname, size=fontsize)
    # do plot
    ax.plot(r[0:i], z[0:i], 'y')
    ax.grid(True)
    # plot a dot if we had a fugitive escape in this interval
    my_escapers.get_new_escapers(t[i])    
    if my_escapers.new_escapers > 0:
        ax.plot(r[i], z[i], marker='o', markersize=2*np.log(my_escapers.new_escapers), markerfacecolor='r', markeredgecolor='r')
    # set the text for the bilions of years ago
    time_in_the_past =  np.max(t) - t[i]
    ax.text(xpos,   ypos,"{:.4f} GYA ".format(time_in_the_past), size=fontsize)
    ax.text(xpos,.8*ypos,"N esc: {:d}".format(my_escapers.total_escapers), size=fontsize)
    


start_time = time.time()
inteval_delay = 100 # miliseconds
print("Creating the animation")
anim = animation.FuncAnimation(fig, animate, frames=my_frames, interval=inteval_delay)
anim.save(outname,dpi=80, writer='imagemagick')
print("--- %s seconds ---" % (time.time() - start_time))
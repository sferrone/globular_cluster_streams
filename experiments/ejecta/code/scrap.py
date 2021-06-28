import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import time 
import h5py
# import astropy.coordinates as coord
import os 
import re

# get backward orbit
GCname='Pal5'
galactic_potential = 'PII'
cluster_potential = 'Plummer'
path = '../../simulations/outputDATA/'+galactic_potential+"/"+cluster_potential+"/"
path_backward = path+"backwardorbits/"+"orbit"+GCname+".dat"
path_streams = path+"streams/"+GCname+".h5"

outname = "../outputs/"+GCname+"_"+galactic_potential+"_"+cluster_potential+"_fugitives.gif"
# the current stream positions 
streams = h5py.File(path_streams, 'r')
# the orbit of the body in time
t, x,y,z,vx,vy,vz = np.loadtxt(path_backward, unpack=True)
# reverse each to make sure we're going in the correct order

# make time go from oldest to now
t = np.flip(t)  # now in hundred of millions of years ago
x = np.flip(x)  
y = np.flip(y)
n_skip = 500
n_frames = int(np.floor(len(x)/n_skip))

# # find out when the stars escape
# fugitives = np.where(streams['energy']['Tescape'][:] > 0)[0]
# time_escape_fugitives = streams['energy']['Tescape'][fugitives] 
# n_fugitives = sum(time_escape_fugitives < t[0])

# First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots(1,1)
fontsize = 15
# set axes limits
the_max = 1.1*np.max([np.max(x),np.max(y)])
the_min = 1.1*np.min([np.min(x),np.min(y)])
ax.set_xlim(the_min, the_max)
ax.set_ylim(the_min ,the_max)
# get coordiantes for witing text
xpos = .675*(the_max-the_min) + the_min
ypos = .925*(the_max-the_min) + the_min
# make graph nice
ax.set_xlabel("Galactic X (kpc)")
ax.set_ylabel("Galactic Y (kpc)")
ax.set_aspect('equal')
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    ax.clear()
    ax.set_xlim(the_min, the_max)
    ax.set_ylim(the_min ,the_max)
    ax.set_xlabel("Galactic X (kpc)", size=fontsize)
    ax.set_ylabel("Galactic Y (kpc)", size=fontsize)    
    ax.set_title(GCname, size=fontsize)
    ax.plot(x[0:n_skip*i], y[0:n_skip*i])
    time_in_the_past = t[n_skip*i]/10
    ax.text(xpos,ypos,"{:.4f} GYA".format(time_in_the_past), size=fontsize)
    return line,


start_time = time.time()
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=n_frames, interval=200, blit=True)
anim.save(outname,dpi=80, writer='imagemagick')
print("--- %s seconds ---" % (time.time() - start_time))



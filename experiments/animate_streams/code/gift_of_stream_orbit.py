# import the relevant libraries
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import animation
from matplotlib.animation import FuncAnimation
from scipy.io import FortranFile
import astropy.coordinates as coord
import astropy.units as u
import time 
import sys
import os 
sys.path.append('../../simulations/code4movies')
import inputMW
# update the style 
plt.rcParams['figure.facecolor'] = 'black'
plt.rcParams['axes.facecolor'] = 'black'
plt.rcParams['axes.grid']=True
plt.rcParams['axes.labelcolor']='w'
plt.rcParams['axes.titlecolor']='w'
plt.rcParams['axes.labelsize']=20
plt.rcParams['grid.alpha'] = .2
plt.rcParams['lines.markerfacecolor']='w'
plt.rcParams['lines.markeredgecolor']='w'
plt.rcParams['scatter.edgecolors']='w'
plt.rcParams['ytick.color']='w'
plt.rcParams['xtick.color']='w'
plt.rcParams['text.color']='w'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['computer modern']
plt.rcParams['font.size']= 15



class Streams:
    """ 
    Store and open the information associated with the stellar streams. 
    We need the name of the globular cluster and the basename of the file and the time step. 
    The time step is assumed to be in units of millions of years 
    """
    def __init__(self,
        GCName,             # i.e. Pal5
        basename            =   "output",   # base name of each bin file
        increment           =   1,          # time step in 10s of millions of years
        current_time        =   0,          # first time step in the sequence
        max_time            =   500,        # this is actually today,
        adjcent_COM_points  =   700         # the number of time stamps to use for the center of mass orbit
        ):
        self.GCName             = GCName
        self.basename           = basename
        self.increment          = increment
        self.current_time       = current_time
        self.max_time           = max_time
        self.gc_frame           = inputMW.ref_frame()
        self.adjcent_COM_points = adjcent_COM_points
        self.folder_path        = \
            "../../simulations/outputDATA4movies/"+GCName+"/"
        self.orbit_file         = \
            "../../simulations/outputDATA/PII/Plummer/backwardorbits/orbit"+GCName+".dat"
        self.set_file_name()
        self.set_center_of_mass_file()
    def set_file_name(self):
        """ Test to see if the stream file exists """ 
        base_time   = str(int(self.current_time)).zfill(3)  # ensure format [102, 213,002]
        fname       = self.folder_path+self.basename+base_time+".bin"
        if not os.path.isfile(fname):
            print("ERROR {} not found".format(fname))
        else:
            self.fname = fname
    def set_center_of_mass_file(self):
        """ Get the file containing the orbit of the center of mass """
        if not os.path.isfile(self.orbit_file):
            print("ERROR!!!: {} not found. \nNo orbital path will be included ".format(
                self.orbit_file))
            self._isCenterOfMassOrbit = False
        else:
            self._isCenterOfMassOrbit = True
            _t_com,_x,_y,_z,_vx,_vy,_vz = \
                np.loadtxt(self.orbit_file, unpack=True) # fancy line continuation
            self._x = np.flip(_x)
            self._y = np.flip(_y)
            self._z = np.flip(_z)
            self._vx = np.flip(_vx)
            self._vy = np.flip(_vy)
            self._vz = np.flip(_vz)

            self._t_com             = _t_com*10 # now in 10 Myrs
            self._n_time_stamps = self._t_com.shape[0]
    def get_COM_trajectory(self):
        """ if the file exists, then get return coordinates near the COM """
        if self._isCenterOfMassOrbit:
            # get the current center of mass position
            com_index=np.argmin(np.abs(self._t_com - self.current_time))
            down_dex = com_index-self.adjcent_COM_points
            up_dex   = com_index+self.adjcent_COM_points
            if down_dex < 0:
                down_dex = 0
            if up_dex > self._n_time_stamps:
                up_dex = self._n_time_stamps
            # extract the COM coordinates
            x  = self._x[down_dex:up_dex]
            y  = self._y[down_dex:up_dex]
            z  = self._z[down_dex:up_dex]
            vx = self._vx[down_dex:up_dex]
            vy = self._vy[down_dex:up_dex]
            vz = self._vz[down_dex:up_dex]
            # store the coordinates of the trajectory in our object
            self.COM_trajectory_Galactocentric = coord.SkyCoord(
                x=x*u.kpc, y=y*u.kpc, z=z*u.kpc, 
                v_x=vx*10*u.km/u.s, v_y=vy*10*u.km/u.s, v_z=vz*10*u.km/u.s, 
                frame=self.gc_frame)
            # store the center of mass coordinates
            self.COM_x = self._x[com_index]
            self.COM_y = self._y[com_index]
            self.COM_z = self._z[com_index]
            self.COM_r = np.sqrt(self._x[com_index]**2 + self._y[com_index]**2)
            self.COM_vx = self._vx[com_index]
            self.COM_vy = self._vy[com_index]
            self.COM_vz = self._vz[com_index]
        else:
            self.COM_trajectory_Galactocentric = None
            self.COM_x = None
            self.COM_y = None
            self.COM_r = None
            self.COM_z = None
            self.COM_vx = None
            self.COM_vy = None
            self.COM_vz = None
    def get_stellar_positions(self):
        # import the positions from the current file
        f = FortranFile(self.fname, 'r')
        x = f.read_reals(dtype='float32')
        y = f.read_reals(dtype='float32')
        z = f.read_reals(dtype='float32')
        vx = f.read_reals(dtype='float32')
        vy = f.read_reals(dtype='float32')
        vz = f.read_reals(dtype='float32')
        f.close() 
        self.Galactocentric = coord.SkyCoord(
            x=x * u.kpc, y=y * u.kpc, z=z * u.kpc, 
            v_x=vx*10*u.km/u.s, v_y=vy*10*u.km/u.s, v_z=vz*10*u.km/u.s, 
            frame=self.gc_frame)
    def increment_system(self):
        # set the new file name, and then get the next set of positions
        self.current_time = self.current_time + self.increment
        self.set_file_name()
        self.get_stellar_positions()
        self.get_COM_trajectory()
    def initialize_plot(self):
        """ start the plot for our system """
        fig, ax = plt.subplots(1,1)
        plt.tight_layout()
        fig.subplots_adjust(bottom=0.2)
        ax.set_aspect('equal')
        self._fig = fig
        self._ax = ax 
    def set_axis_properties(self,
        x_min=-20,x_max=20,y_min=-20,y_max=20,xlabel="x", ylabel="y"):
        """ set the generic plot properties here """
        font_size = 15
        # assumes units of kpc
        self._ax.set_xlim(x_min, x_max)
        self._ax.set_ylim(y_min ,y_max)
        self._ax.set_xlabel(xlabel, size=font_size)
        self._ax.set_ylabel(ylabel, size=font_size)
        self._ax.tick_params(axis='x', labelrotation = 30)
        self._ax.grid(True)


# define the animate function 
def animate(i, my_stream, fontsize, xpos, ypos, tmax=500):
    my_stream._ax.clear()
    my_stream.increment_system()
    # plot the trajectory and center of mass
    if my_stream._isCenterOfMassOrbit:
        my_stream._ax.plot(my_stream.COM_trajectory_Galactocentric.x, 
        my_stream.COM_trajectory_Galactocentric.y,
        color='r', zorder=1, alpha=.5)
        my_stream._ax.scatter(my_stream.COM_x, my_stream.COM_y, color='r', zorder=2, alpha=.8)
    my_stream._ax.scatter(my_stream.Galactocentric.x,my_stream.Galactocentric.y, 
    s=.05, alpha=.8, zorder=3)
    my_stream.set_axis_properties(xlabel="x (kpc)", ylabel="y (kpc)")
    time_in_the_past = (tmax - i) / 100 # put 10s of millions of years into Gyrs
    my_stream._ax.text(xpos,   ypos,"{:.4f} GYA ".format(time_in_the_past), size=fontsize)


def main():
    """ set the parameters 
    run the animation """

    # in the main you have to start running this stuff
    # create the object
    GCName              = "Pal5"
    my_stream           = Streams(GCName)
    galactic_potential  = "PII"
    GCpotential         = "Plummer"
    coordainte_system   = "galactocentric"
    units               = "xy"
    outname             = \
        "../outputs/"+coordainte_system+"_"+galactic_potential+\
            "_"+GCpotential+"_"+units+"_"+GCName+".gif"
    # establish the proper number of iterations
    start_time = 1
    end_time = 30 # in 10s of millions of years 
    interval = 1
    my_frames = range(start_time,end_time,interval)
    # set properties
    delay = 100 # miliseconds
    fontsize = 15
    x_min = -20
    x_max =  20
    y_min = -20
    y_max =  20    
    xpos = .675*(x_max-x_min) + x_min
    ypos = .925*(y_max-y_min) + y_min

    my_stream.initialize_plot()
    
    start_time = time.time()
    anim = animation.FuncAnimation(my_stream._fig,animate,frames=my_frames,interval=delay,
        fargs=(my_stream, fontsize, xpos, ypos))

    anim.save(outname,dpi=150, writer='imagemagick')
    print("--- %s seconds ---" % (time.time() - start_time))
if __name__=="__main__":
    main()
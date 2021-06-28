import numpy as np 
import matplotlib.pyplot as plt 
import h5py
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u 
from matplotlib.colors import LogNorm
import sys
from PIL import Image, ImageOps, ImageDraw, ImageFont

sys.path.append('../../functions')
import global_functions as my_gf 

def get_stellar_surface_density_galactic_coordinates(galactic_coordinates, n_bins = 500):
    """ get the surface density for the plot """
    # shift the zero longitude [-180 to 180]
    lons = galactic_coordinates.l.value
    lons[lons > 180] = lons[lons > 180] - 360
    lats = galactic_coordinates.b.value
    lon_limits = [-180,180]
    lat_limits = [-91,91]
    # get the counts per bin
    rangev=[lon_limits,lat_limits]
    count, xedges, yedges =np.histogram2d(lons, lats, bins=n_bins,range=rangev)
    count = np.rot90(count)
    count = np.flipud(count)
    # create the meshgrid 
    xv,yv = np.meshgrid(np.linspace(-180,180,n_bins),np.linspace(-92,92,n_bins))
    # only take bins with data
    x_flat = xv[count > 0]
    y_flat = yv[count > 0]
    z_flat = count[count > 0]
    # ensure that points of higher density are shown on top 
    idx = z_flat.argsort()
    x_flat, y_flat, z_flat = x_flat[idx], y_flat[idx], z_flat[idx]
    # put into counts per stradian
    z_flat *= n_bins/(4*np.pi)
    return x_flat, y_flat, z_flat

def get_xy_stellar_surface_density_galactocentric_coordinates(galactocentric_coords, n_bins = 500):
    """ get the stellar surface density in units of counts per kpc^2 """
    # store the values more conviently
    my_x = galactocentric_coords.x.value
    my_y = galactocentric_coords.y.value
    # we want out limits in the xy plane to be square
    upper_limit = np.ceil(np.max(abs(np.concatenate([my_y, my_x]))))
    limits = [-upper_limit, upper_limit]
    rangev=[limits,limits]
    count, xedges, yedges =np.histogram2d(my_x, my_y, bins=n_bins,range=rangev)
    count = np.rot90(count)
    count = np.flipud(count)
    # create the meshgrid
    spatial_array = np.linspace(limits[0],limits[1],n_bins)
    xv,yv = np.meshgrid(spatial_array,spatial_array)
    # only take bins with data
    x_flat = xv[count > 0]
    y_flat = yv[count > 0]
    z_flat = count[count > 0]
    # ensure that points of higher density are shown on top 
    idx = z_flat.argsort()
    x, y, planar_stellar_surface_density = x_flat[idx], y_flat[idx], z_flat[idx]
    # place into the proper units
    total_area = (spatial_array[1]-spatial_array[0])*(spatial_array[1]-spatial_array[0])
    planar_stellar_surface_density *= n_bins / total_area
    return x, y, planar_stellar_surface_density
    
def get_rz_stellar_surface_density_galactocentric_coordinates(galactocentric_coords, n_bins = 500):
    """ get the stellar surface density in units of counts per kpc^2 """
    # store the values more conviently
    my_x = galactocentric_coords.x.value
    my_y = galactocentric_coords.y.value
    my_z = galactocentric_coords.z.value
    my_r = np.sqrt( my_x**2 + my_y**2 ) 
    # we want out limits in the xy plane to be square
    upper_limit = np.ceil(np.max(my_r))
    r_limits = [0, upper_limit]
    z_limit = np.ceil(np.max(abs(my_z)))
    rangev=[r_limits,[-z_limit,z_limit]]
    count, xedges, yedges =np.histogram2d(my_r, my_z, \
        bins=n_bins,range=rangev)
    count = np.rot90(count)
    count = np.flipud(count)
    # create the meshgrid
    z_space = np.linspace(-z_limit,z_limit,n_bins)
    r_space = np.linspace(r_limits[0],r_limits[1],n_bins)
    xv,yv = np.meshgrid(r_space,z_space)
    # only take bins with data
    r_flat = xv[count > 0]
    z_flat = yv[count > 0]
    count_flat = count[count > 0]
    # ensure that points of higher density are shown on top 
    idx = count_flat.argsort()
    r, z, cylindrical_surface_density = r_flat[idx], z_flat[idx], count_flat[idx]
    # size per bin
    bin_size = (z_space[1]-z_space[0])*(r_space[1]-r_space[0]) / n_bins
    cylindrical_surface_density /= bin_size
    return r, z, cylindrical_surface_density

def plot_xy_rz_stellar_surface_density(cluster_name, x,y,planar_stellar_surface_density, r, z, cylindrical_surface_density):
    """ do the plot like your life depends on it """
    # get some preliminary information
    cbar_name   = r'$\Sigma^{{*}}$ (counts / $kpc^2$)'
    upper_limit = np.ceil(np.max(abs(np.concatenate([x, y]))))
    limits      = [-upper_limit, upper_limit]
    out_name    = cluster_name+"_xy_rz.png"
    my_fontsize = 20
    dot_size    = 1
    # star the plot 
    fig, ax = plt.subplots(1,2)
    fig.set_figheight(5)
    fig.set_figwidth(11)    
    fig.subplots_adjust(wspace=.05, hspace=0)
    # plot rz
    im2 = ax[0].scatter(r,z,c=cylindrical_surface_density,\
        s=dot_size, cmap="rainbow",norm=LogNorm())
    # add the colorbar
    cbar = fig.colorbar(im2, ax=ax,fraction=0.046, pad=0.1)
    cbar.ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
    cbar.ax.get_yaxis().labelpad = 20
    # plot xy
    im = ax[1].scatter(x, y, c=planar_stellar_surface_density, \
        s=dot_size, cmap="rainbow", norm=LogNorm(vmin=cbar.vmin, vmax=cbar.vmax))
    # take care of the axis limits
    z_max       = np.ceil(np.max(abs(z)))
    r_max       = np.ceil(np.max(r))
    # make sure that the axis are square and cover the whole data range
    if 2*z_max > r_max:
        z_limits    = [-z_max, z_max]
        r_limits    = [0, 2*z_max]
    elif r_max >= 2*z_max:
        z_max       = np.ceil(r_max/2)
        r_max       = 2*z_max
        z_limits    = [-z_max, z_max]
        r_limits    = [0, r_max]
    # set plot properties
    ax[0].grid(True, alpha=.2)
    ax[0].set_xlim(r_limits)
    ax[0].set_ylim(z_limits)
    ax[0].set_aspect('equal', 'box')
    # plot planar
    ax[1].grid(True, alpha=.2)
    ax[1].set_xlim(limits)
    ax[1].set_ylim(limits)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].set_aspect('equal', 'box')
    # set the label information
    ax[0].set_xlabel("R (kpc)", size=.75*my_fontsize)
    ax[0].set_ylabel("Z (kpc)", size=.75*my_fontsize)
    ax[1].set_xlabel("X (kpc)", size=.75*my_fontsize)
    ax[1].set_ylabel("Y (kpc)", size=.75*my_fontsize)
    # add the position of the sun 
    ax[0].plot(8,0, '*',  color='w', markersize=12, markeredgecolor='k')
    ax[1].plot(-8,0, '*', color='w', markersize=12, markeredgecolor='k')
    # plot the center of mass
    glob_center,found_names = my_gf.get_COM_coordinates_RA([cluster_name])
    plot_COM_position(ax[0],glob_center, coordinate_system="rz")
    plot_COM_position(ax[1],glob_center, coordinate_system="xy")
    fig.savefig("../density_plots/"+out_name, dpi=150,bbox_inches='tight')

def plot_stellar_surface_density_galactic_coordinates(\
    cluster_name, lons, lats, stellar_surface_density):
    """ do the plot and make it nice """
    # limits for the axis
    cbar_name = r'$\Sigma^{{*}}$ (counts / steradian)'
    xlims = [-180, 180]
    ylims = [-90, 90]
    # other plot properties
    my_fontsize = 20
    dot_size    = 1
    fig, ax     = plt.subplots(1,1)
    fig.set_figheight(5)
    fig.set_figwidth(11)
    # plot the data
    im = ax.scatter(lons, lats, c=stellar_surface_density, s=dot_size, norm=LogNorm(), cmap="rainbow")
    # add the colorbar
    cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
    cbar.ax.get_yaxis().labelpad = 20
    # make the axes nice
    ax.grid(True, alpha=.2)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_aspect('equal')
    # make consistent with ibata 2021
    ax.invert_xaxis()
    ax.set_xticks(np.arange(-180,210,30))
    ax.set_yticks(np.arange(-90,120,30))
    # add the plot information 
    ax.set_xlabel("Longitude",   size=my_fontsize)
    ax.set_ylabel("Latitude",    size=my_fontsize)
    out_name = cluster_name+"_galactic_coords.png"
    fig.savefig("../density_plots/"+out_name, dpi=150,bbox_inches='tight')
    return cbar

def plot_backward_orbit(cluster_name,tb,xb,yb,zb):
    """ we want to color this in time """
    # get the radius
    rb = np.sqrt(xb**2 + yb**2)
    # get the RA and DEC object
    print(cluster_name)
    glob_center,found_names = my_gf.get_COM_coordinates_RA([cluster_name])
    print("found_names", found_names)
    # sort to ensure that the most recent is on top
    idx = tb.argsort()
    idx = idx[::-1]
    tb, xb, yb, zb, rb = tb[idx], xb[idx], yb[idx], zb[idx], rb[idx]
    # some plot stuff
    cbar_name = "Look back time (GYA)"
    upper_limit = np.ceil(np.max(abs(np.concatenate([xb, yb]))))
    limits      = [-upper_limit, upper_limit]
    out_name    = "lookback_orbit_xy_rz{:s}.png".format(cluster_name)
    my_dot_size = 1
    my_fontsize = 20
    # start the plot 
    fig, ax = plt.subplots(1,2)
    fig.set_figheight(5)
    fig.set_figwidth(11)    
    fig.subplots_adjust(wspace=.05, hspace=0)
    # plot rz
    im = ax[0].scatter(rb,zb,s=my_dot_size,c=tb,cmap='jet_r',vmin=0.,vmax=5.)
    # add the colorbar
    cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.1)
    cbar.ax.set_ylabel(cbar_name, rotation=270, size=.75*my_fontsize)
    cbar.ax.get_yaxis().labelpad = 20
    # plot xy
    im2 = ax[1].scatter(xb,yb,s=my_dot_size,c=tb,cmap='jet_r',vmin=cbar.vmin,vmax=cbar.vmax)
    # take care of the axis limits
    z_max       = np.ceil(np.max(abs(zb)))
    r_max       = np.ceil(np.max(rb))
    # make sure that the axis are square and cover the whole data range
    if 2*z_max > r_max:
        z_limits    = [-z_max, z_max]
        r_limits    = [0, 2*z_max]
    elif r_max <= 2*z_max:
        z_max       = np.ceil(r_max/2)
        r_max       = 2*z_max
        z_limits    = [-z_max, z_max]
        z_limits    = [0, r_max]
    # set plot properties
    ax[0].grid(True, alpha=.2)
    ax[0].set_xlim(r_limits)
    ax[0].set_ylim(z_limits)
    ax[0].set_aspect('equal', 'box')
    # plot planar
    ax[1].grid(True, alpha=.2)
    ax[1].set_xlim(limits)
    ax[1].set_ylim(limits)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].set_aspect('equal', 'box')
    # set the label information
    ax[0].set_xlabel("R (kpc)", size=.75*my_fontsize)
    ax[0].set_ylabel("Z (kpc)", size=.75*my_fontsize)
    ax[1].set_xlabel("X (kpc)", size=.75*my_fontsize)
    ax[1].set_ylabel("Y (kpc)", size=.75*my_fontsize)
    # add the position of the sun 
    ax[0].plot(8,0, '*',  color='w', markersize=12, markeredgecolor='k')
    ax[1].plot(-8,0, '*', color='w', markersize=12, markeredgecolor='k')
    # plot the center of mass
    plot_COM_position(ax[0], glob_center, coordinate_system="rz")
    plot_COM_position(ax[1], glob_center, coordinate_system="xy")

    fig.savefig("../density_plots/"+out_name, dpi=150,bbox_inches='tight')
    
def plot_COM_position(ax, glob_center, coordinate_system="galactic_coords"):
    """ plot the center of mass """
    if coordinate_system=="galactic_coords":
        galactic_coordinates = glob_center.transform_to(coord.Galactic)
        my_x = galactic_coordinates.lon.value
        my_y = galactic_coordinates.lat.value
        # make between [-180,180]
        my_x[my_x > 180] = my_x[my_x > 180] - 360
    elif coordinate_system=="xy":
        galactocentric = glob_center.transform_to(coord.Galactocentric)
        my_x = galactocentric.x.value
        my_y = galactocentric.y.value
    elif coordinate_system=="rz":
        galactocentric = glob_center.transform_to(coord.Galactocentric)
        x = galactocentric.x.value
        y = galactocentric.y.value
        z = galactocentric.z.value
        r = np.sqrt(x**2 + y**2)
        my_x = r 
        my_y = z
    print("(x,y) = ({},{})".format(my_x, my_y))
    ax.scatter(my_x,my_y,c='w',s=20,edgecolors="k")

def concatenate_images(cluster_name):
    """ concatenate the produced images into one larger figure """
    # set path infromation  
    base_path = '../density_plots/'
    galactic_coords = cluster_name+"_galactic_coords.png"
    xy_rz_coords    = cluster_name+"_xy_rz.png"
    look_back       = "lookback_orbit_xy_rz{:s}.png".format(cluster_name)
    image_paths = [base_path+galactic_coords, base_path+xy_rz_coords, base_path+look_back ]
    # open the images 
    imgs = [Image.open(i) for i in image_paths]
    # pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
    min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
    imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )

    # for a vertical stacking it is simple: use vstack
    imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
    imgs_comb = Image.fromarray(imgs_comb)
    imgs_comb.save('../density_plots/trifecta_'+cluster_name+".png")


def main():
    # Import some data
    potential_milk  = "PII"
    simulation_path = "../../simulations/outputDATA/"+potential_milk+"/streams/"
    backward_o_path = "../../simulations/outputData/"+potential_milk+"/backwardorbits/orbit"
    cluster_name = "Pal5"
    # extract data
    x,y,z,vx,vy,vz=my_gf.get_stream(simulation_path,cluster_name,5)
    tb,xb,yb,zb,vxb,vyb,vzb=my_gf.get_backward_orbit(backward_o_path+cluster_name)
    # put into astropy structure
    galactocentric_coords = coord.Galactocentric(
        x = x * u.kpc,
        y = y * u.kpc,
        z = z * u.kpc,
        v_x = vx * 10 * u.km / u.s,
        v_y = vy * 10 * u.km / u.s,
        v_z = vz * 10 * u.km / u.s,
    )
    # transform to galactic coordinates
    galactic_coordinates = galactocentric_coords.transform_to(coord.Galactic)
    # get coordinates  
    lons,lats,SSD   = get_stellar_surface_density_galactic_coordinates(galactic_coordinates)
    x,y,planar      = get_xy_stellar_surface_density_galactocentric_coordinates(galactocentric_coords, n_bins = 500)
    r,z,cylin       = get_rz_stellar_surface_density_galactocentric_coordinates(galactocentric_coords, n_bins = 500)
    # plot in galactic coordinates
    plot_stellar_surface_density_galactic_coordinates(cluster_name, lons,lats,SSD)
    plot_xy_rz_stellar_surface_density(cluster_name, x,y,planar, r,z,cylin)
    plot_backward_orbit(cluster_name, tb, xb, yb, zb)
    concatenate_images(cluster_name)
if __name__=="__main__":
    main()


# /mod4gaia/GCsTT/outputDATA/PII
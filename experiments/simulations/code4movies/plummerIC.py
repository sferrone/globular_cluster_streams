"""
 GENERATE POSITIONS AND VELOCITIES OF A PLUMMER SPHERE
 by Montacarlo sampling
"""

import numpy as np
import matplotlib.pyplot as plt
import random as random
import time
import inputGC
from astropy.io import fits
from astropy.table import Table

"""
Assign parameters: core radius, total mass, fraction of mass inside cutoff radius, total number of particles
"""

f = open('GClist', 'r')
for line in f :
    GCname = line.strip()
    print(GCname)

    G=1.
    mtot,rcore=inputGC.Plummer(GCname)
    cutoff_massfraction=0.990 # fraction of mass inside the cutoff radius
    ntot=100000 # number of particles to be generated

    mass=[cutoff_massfraction*mtot/ntot]*ntot # mass of a single particle
    cutoff_mass=cutoff_massfraction*mtot
    c=cutoff_massfraction**(2/3)
    cutoff_radius=np.sqrt(rcore*rcore*(c/(1-c)))
    print("cutoff radius "+str(cutoff_radius))

    """
    Define functions
    """
    def plummer(mtot,rcore,r):
        plummer=3*mtot/(4*np.pi*rcore**3)
        plummer=plummer/(1+(r/rcore)**2)**(2.5)
        return plummer

    def df(mtot,rcore,x,vv):
        E=0.5*vv*vv-G*mtot/((x*x+rcore*rcore)**0.5)
        if(E > 0):
            print(E)
            print("error")
        df=(-E)**3.5
        return df

    def escape_velocity(xr,mtot,rcore):

        vesc=2*np.abs(-G*mtot/np.sqrt(xr*xr+rcore*rcore))
        vesc=np.sqrt(vesc)

        return vesc


    """
    Generate positions and velocities
    """

    start=time.time()    

    dr=0.1*rcore
    r=np.arange(0,cutoff_radius,cutoff_radius*dr) # define radii
    plummertot=[plummer(mtot,rcore,r[i]) for i in range(0,len(r))] # plummer density as a function of r
    maxf=np.amax(plummertot*r**2) # calculate maximum of the plummertot*r**2 function


    count=0
    radius=[]
    velocity=[] 
    yextr=[]
    xrej=[]
    yrej=[]

    while count < ntot:
        x=np.random.uniform(0,cutoff_radius,1)
        rho=np.random.uniform(0,maxf,1)
        fx=np.interp(x, r,plummertot*r**2)
        if(rho <= fx):
            count+=1
            radius.append(x)
            yextr.append(rho)
            countv=0
            while countv < 1:
                vesc=escape_velocity(x,mtot,rcore)
                dv=vesc*0.01
                v=np.arange(0,vesc,dv)
                dftot=[df(mtot,rcore,float(x),v[i]) for i in range(0,len(v))]
                maxfv=np.amax(dftot*v**2)
                xv=np.random.uniform(0,vesc,1)
                yv=np.random.uniform(0,maxfv,1)
                fxv=np.interp(xv, v, dftot*v**2)
                if(yv <= fxv):
                    countv+=1
                    velocity.append(xv)
        else:
            xrej.append(x)
            yrej.append(rho)

    theta=np.random.uniform(0,2*np.pi,ntot)
    cosphi=np.random.uniform(-1.,1.,ntot)
    phi=np.arccos(cosphi)

    radius=np.concatenate(radius, axis=0 )
    
    x=np.multiply(np.multiply(radius,np.cos(theta)),np.sin(phi))
    y=np.multiply(np.multiply(radius,np.sin(theta)),np.sin(phi))
    z=np.multiply(radius,np.cos(phi))

    theta=np.random.uniform(0,2*np.pi,ntot)
    cosphi=np.random.uniform(-1.,1.,ntot)
    phi=np.arccos(cosphi)

    velocity=np.concatenate( velocity, axis=0 )
    vx=np.multiply(np.multiply(velocity,np.cos(theta)),np.sin(phi))
    vy=np.multiply(np.multiply(velocity,np.sin(theta)),np.sin(phi))
    vz=np.multiply(velocity,np.cos(phi))

    print(vx)
    print(vy)
    print(vz)
    
    end=time.time()
    print("elapsed time "+str(end-start))


    """
    Save output file
    """

    tbdata=Table([x,y,z,vx,vy,vz],names=('x0','y0','z0','vx0','vy0','vz0'))
    tbdata.write("../plummer_profiles/"+GCname+'Plummer.fits', format='fits', overwrite=True)

    """
    Figures to check that the generation of the positions is correct
    """
    plt.figure()
    plt.plot(r,plummertot*r**2,label='Theoretical curve')
    plt.yscale('log')
    plt.xscale('log')
    plt.axhline(maxf,color='red')

    plt.scatter(xrej,yrej,color='grey',s=2,label='Rejected radii')
    plt.scatter(radius,yextr,color='red',s=5,label='Accepted  radii')
    plt.legend()
    plt.title(GCname, size=15)
    plt.savefig("../plummer_profiles/"+"density_profiles/"+GCname+"profile.png", dpi=150)
    plt.figure()
    plt.yscale('log')
    plt.xscale('log')


    dr=0.02*rcore
    imax=np.int(cutoff_radius/dr)
    icount=[0]*(imax+1)
    for i in range(0,ntot):
        ix=np.int(radius[i]/dr)    
        icount[ix]+=mass[i]

    xp=[]
    yp=[]
    for ix in range(1,imax+1):
        if(icount[ix] > 0):
            icount[ix]=icount[ix]/(4.*np.pi*dr*(ix*dr)**2)
            xp.append(ix*dr)
            yp.append(icount[ix])

    plt.plot(xp,yp,label='Density profile from generated positions')
        
    plt.plot(r,plummertot,color='k',label='Theoretical density profile')
    plt.legend()
    plt.title(GCname, size=15)
    plt.savefig("../plummer_profiles/"+"density_profiles/"+GCname+"density_profile.png", dpi=150)
    #plt.show()
    plt.close('all')

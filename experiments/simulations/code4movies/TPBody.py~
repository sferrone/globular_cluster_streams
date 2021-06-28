import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u
import time
import inputMW
import inputGC
import output
import os
from integrate import *
import h5py

"""                                                                                                        
GALACTIC REFERENCE FRAME
"""
gc_frame=inputMW.ref_frame()

"""
ASSIGN MILKY WAY MODEL
"""

model='PII'
modelshort=model


"""                    
ASSIGN GLOBULAR CLUSTER AND GET ITS GALACTOCENTRIC POSITIONS AND VELOCITIES
"""
n_trials=1 # for monte carlo error simulation

f = open('GClistred', 'r')
for line in f :
    GCname = line.strip()
    print(GCname)
    
    for number in range(0,n_trials):
        ierror=0
        if(number > 0):
            ierror=1

        print('=================')
        if(ierror==0):
            print('SIMULATION WITHOUT ERRORS')
        if(ierror==1):
            print('SIMULATION WITH ERRORS # '+str(number))

        xGC,yGC,zGC,vxGC,vyGC,vzGC=inputGC.GCcoord(GCname,gc_frame,ierror)
        print("current x,y,z,vx,vy,vz")
        print(xGC,yGC,zGC)
        print(vxGC,vyGC,vzGC)


        """
        ASSIGN TOTAL NUMBER OF STEP FOR THE ORBIT INTEGRATION AND INTEGRATION TIME-STEP 
        """

        Nstep=50000
        deltat=0.001
        isave=100 # save a snapshot every isave steps

        """
        INTEGRATE THE ORBITS OF THE GC BARYCENTRE BACKWARD IN TIME, TO GET THE "INITIAL CONDITIONS"
        """

        backward='YES'
        N=1

        output.createdir('../outputDATA/')
        pathorbits='/backwardorbits'
        output.createdir('../outputDATA/'+model)
        output.createdir('../outputDATA/'+model+pathorbits)
        if(ierror==1):
            output.createdir('../outputDATA/'+model+pathorbits+'/err/')
            output.createdir('../outputDATA/'+model+pathorbits+'/err/'+GCname)
            output.createdir('../outputDATA/'+model+pathorbits+'/err/'+GCname+'/'+str(number))
            pathorbits='/backwardorbits/err/'+GCname+'/'+str(number)
        GCparam=[0.,0.]
        MWparam=inputMW.MW(modelshort)
        integrate.orbits(xGC,yGC,zGC,vxGC,vyGC,vzGC,GCparam,MWparam,model,modelshort,GCname,backward,pathorbits,deltat,Nstep,isave,N)

        
        xGCini=integrate.xp
        yGCini=integrate.yp
        zGCini=integrate.zp
        vxGCini=integrate.vxp
        vyGCini=integrate.vyp
        vzGCini=integrate.vzp

        print("initial x,y,z,vx,vy,vz")
        print(xGCini,yGCini,zGCini)
        print(vxGCini,vyGCini,vzGCini)


        """                                                                                                        
        READ INITIAL CONDITIONS OF THE N-BODY GC
        """
        x0,y0,z0,vx0,vy0,vz0=inputGC.Nbody(GCname)


        #then place the GC in the Galaxy, assigning to its baricentre positions and velocities as T Gyr ago 
        x0=x0+xGCini
        y0=y0+yGCini
        z0=z0+zGCini
        vx0=vx0+vxGCini
        vy0=vy0+vyGCini
        vz0=vz0+vzGCini
    
    
        #now integrate the N-body GC forward in time
        start=time.time()

        backward='NO'
        N=len(x0)    

        GCparam=inputGC.Plummer(GCname)
        MWparam=inputMW.MW(modelshort)
    
        integrate.deallocate()
        integrate.orbits(x0,y0,z0,vx0,vy0,vz0,GCparam,MWparam,model,modelshort,GCname,backward,pathorbits,deltat,Nstep,isave,N)

        xf=integrate.xp
        yf=integrate.yp
        zf=integrate.zp
        vxf=integrate.vxp
        vyf=integrate.vyp
        vzf=integrate.vzp
        tesc=integrate.tesc
        phiMW=integrate.phimw
        phiGC=integrate.phigc

        end=time.time()
        print("elapsed time "+str(end-start))

        c1 = coord.SkyCoord(x=xf * u.kpc, y=yf * u.kpc,z=zf * u.kpc, v_x=vxf*10 * u.km/u.s, v_y=vyf*10 * u.km/u.s, v_z=vzf*10 * u.km/u.s, frame=gc_frame)
        gc2=c1.transform_to(coord.ICRS) 
    
        ra=gc2.ra.value
        dec=gc2.dec.value
        D=gc2.distance.value
        pm_ra_cosdec=gc2.pm_ra_cosdec.value
        pm_dec=gc2.pm_dec.value
        RV=gc2.radial_velocity.value

        gc3=gc2.transform_to('galactic') 
        ll=gc3.l.value
        l=ll
        l[ll>180]=l[ll>180]-360.
        b=gc3.b.value
        pm_l_cosb=gc3.pm_l_cosb.value
        pm_b=gc3.pm_b.value


        """                                                                                                        
        SAVE OUTPUT FILE
        """
    

        if(ierror==0):
            pathorbits='/streams/'
            output.createdir('../outputDATA/'+model+pathorbits)
            hf = h5py.File('../outputDATA/'+model+pathorbits+GCname+'.h5', 'w')
        if(ierror==1):
            output.createdir('../outputDATA/'+model+pathorbits+GCname)
            output.createdir('../outputDATA/'+model+pathorbits+GCname+'/err/')
            pathorbits='/streams/'+GCname+'/err/'
            hf = h5py.File('../outputDATA/'+model+pathorbits+GCname+str(number)+'.h5', 'w')
        EQ = hf.create_group('equatorial')
        EQ.create_dataset('RA',data=ra)
        EQ.create_dataset('DEC',data=dec)
        EQ.create_dataset('D',data=D)
        EQ.create_dataset('PMRA_COSDEC',data=pm_ra_cosdec)
        EQ.create_dataset('PMDEC',data=pm_dec)
        EQ.create_dataset('RV',data=RV)
        GAL = hf.create_group('galactic')
        GAL.create_dataset('LONG',data=l)
        GAL.create_dataset('LAT',data=b)
        GAL.create_dataset('PML_COSB',data=pm_l_cosb)
        GAL.create_dataset('PMB',data=pm_b)
        GALCEN = hf.create_group('galactocentric')
        GALCEN.create_dataset('X',data=xf)
        GALCEN.create_dataset('Y',data=yf)
        GALCEN.create_dataset('Z',data=zf)
        GALCEN.create_dataset('VX',data=vxf)
        GALCEN.create_dataset('VY',data=vyf)
        GALCEN.create_dataset('VZ',data=vzf)
        ENERGY = hf.create_group('energy')
        ENERGY.create_dataset('Tescape',data=tesc)
        ENERGY.create_dataset('phiMW',data=phiMW)
        ENERGY.create_dataset('phiGC',data=phiGC)

        hf.close()

        integrate.deallocate()
        print("END")    




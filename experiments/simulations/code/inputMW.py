import numpy as np
import astropy.coordinates as coord
import astropy.units as u

"""    
Reference frame: values
"""    
def ref_frame():

    vLSR= [0., 240.,0] * (u.km / u.s) # from Reid et al 2014   
    v_sun = [11.1, 12.24, 7.25] * (u.km / u.s)  # [vx, vy, vz] from Schonrich (2012)                                                    
    gc_frame = coord.Galactocentric(galcen_distance = 8340*u.pc,galcen_v_sun=vLSR+v_sun,z_sun=27*u.pc) # z-sun from Chen et al. 2001, galcen_distance from Reid et al 2014  

    return gc_frame


"""    
Galaxy model parameters
"""    
def MW(modelshort):

    if(modelshort=="AS"):
        Md_bulge=606.
        Md_halo=4615.
        Md_thin=3690.
        Md_thick=0.
        Md_bar=0
        zd_bulge=0.25
        rd_bulge=0.
        axhalo = 12.
        bxhalo = 0.
        cxhalo = 0.
        thetahalo = 0.
        zd_thin=0.25
        rd_thin=5.3178
        zd_thick=0.
        rd_thick=0.
        axbar=0.
        bxbar=0.
        cxbar=0.

    if(modelshort=="PI"):
        Md_bulge=460.
        Md_halo=6000.
        Md_thin=1700.
        Md_thick=1700.
        Md_bar=0.
        zd_bulge=0.3
        rd_bulge=0.
        axhalo = 14.
        bxhalo = 0.
        cxhalo = 0.
        thetahalo = 0.
        zd_thin=0.25
        rd_thin=5.3
        zd_thick=0.8
        rd_thick=2.6
        axbar=0.
        bxbar=0.
        cxbar=0.


    if(modelshort=="PII"):
        Md_bulge=0.
        Md_halo=9000.
        Md_thin=1600.
        Md_thick=1700.
        Md_bar=0.
        zd_bulge=0.
        rd_bulge=0.
        axhalo = 14.
        bxhalo = 0.
        cxhalo = 0.
        thetahalo = 0.
        zd_thin=0.25
        rd_thin=4.8
        zd_thick=0.8
        rd_thick=2.0
        axbar=0.
        bxbar=0.
        cxbar=0.

    if(modelshort=="triaxial"):
        Md_bulge=0.
        Md_halo=14000.
        Md_thin=1600.*1.4
        Md_thick=1700.*1.4
        Md_bar=0.
        zd_bulge=0.
        rd_bulge=0.
        axhalo = 10.
        bxhalo = 10.
        cxhalo = 10.
        thetahalo = 0.
        zd_thin=0.25
        rd_thin=4.8
        zd_thick=0.8
        rd_thick=2.0
        axbar=0.
        bxbar=0.
        cxbar=0.


    if(modelshort=="PII_0.2_FAST"):
        Md_bulge=0.
        Md_halo=9000.
        Md_thin=1600.
        Md_thick=1700.
        Md_bar=0.2*(Md_thin+Md_thick)
        Md_thin=0.8*Md_thin
        Md_thick=0.8*Md_thick
        zd_bulge=0.
        rd_bulge=0.
        zd_halo=14.
        axhalo = 14.
        bxhalo = 0.
        cxhalo = 0.
        thetahalo = 0.        
        zd_thin=0.25
        rd_thin=4.8
        zd_thick=0.8
        rd_thick=2.0
        axbar=2.
        bxbar=1.
        cxbar=0.5

    if(modelshort=="PII_0.3_FAST"):
        Md_bulge=0.
        Md_halo=9000.
        Md_thin=1600.
        Md_thick=1700.
        Md_bar=0.3*(Md_thin+Md_thick)
        Md_thin=0.7*Md_thin
        Md_thick=0.7*Md_thick
        zd_bulge=0.
        rd_bulge=0.
        axhalo = 14.
        bxhalo = 0.
        cxhalo = 0.
        thetahalo = 0.
        zd_thin=0.25
        rd_thin=4.8
        zd_thick=0.8
        rd_thick=2.0
        axbar=2.
        bxbar=1.
        cxbar=0.5

    if(modelshort=="PII_0.2_SLOW"):
        Md_bulge=0.
        Md_halo=9000.
        Md_thin=1600.
        Md_thick=1700.
        Md_bar=0.2*(Md_thin+Md_thick)
        Md_thin=0.8*Md_thin
        Md_thick=0.8*Md_thick
        zd_bulge=0.
        rd_bulge=0.
        axhalo = 14.
        bxhalo = 0.
        cxhalo = 0.
        thetahalo = 0.
        zd_thin=0.25
        rd_thin=4.8
        zd_thick=0.8
        rd_thick=2.0
        axbar=4.
        bxbar=1.
        cxbar=0.5


    if(modelshort=="PII_0.3_SLOW"):
        Md_bulge=0.
        Md_halo=9000.
        Md_thin=1600.
        Md_thick=1700.
        Md_bar=0.3*(Md_thin+Md_thick)
        Md_thin=0.7*Md_thin
        Md_thick=0.7*Md_thick
        zd_bulge=0.
        rd_bulge=0.
        axhalo = 14.
        bxhalo = 0.
        cxhalo = 0.
        thetahalo = 0.
        zd_thin=0.25
        rd_thin=4.8
        zd_thick=0.8
        rd_thick=2.0
        axbar=4.
        bxbar=1.
        cxbar=0.5


    MWparam=[Md_bulge,Md_halo,Md_thin,Md_thick,Md_bar,zd_bulge,rd_bulge,axhalo,bxhalo,cxhalo,thetahalo,zd_thin,rd_thin,zd_thick,rd_thick,axbar,bxbar,cxbar]

    return MWparam


import numpy as np
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import fits

""" GC parameters (Plummer sphere) """
fname_baumgardt = "../cluster_data/baumgardt_structural_data.fits"
fname_vasiliev = "../cluster_data/asu_with_baumgardt_names.fit"

def Plummer(GCname):
    fits_baumgardt = fits.open(fname_baumgardt)
    cluster_index = fits_baumgardt[1].data['ID']==GCname
    if(sum(cluster_index)==0):
        print("ERROR: Cluster", GCname,"not found in", fname_baumgardt)
    elif (sum(cluster_index) > 1):
        print("ERROR: Cluster", GCname,"Found multiple times in", fname_baumgardt)
    else:
        half_mass_radius = fits_baumgardt[1].data['rh_m'][cluster_index][0]
        a = ((1/2)**(-2/3) -1)**(1/2) * half_mass_radius / 1000 # characteristic radius in kpc
        Mtot = fits_baumgardt[1].data['Mass'][cluster_index][0] * 10**-7 # total mass                                
    MGCparam=[Mtot,a]
    fits_baumgardt.close()
    return MGCparam

def allPlummer(allGCname):

    allGCparam=[]
    fits_baumgardt = fits.open(fname_baumgardt)
    for i in range(0,len(allGCname)):
        GCname = allGCname[i]
        cluster_index = fits_baumgardt[1].data['ID']==GCname
        if(sum(cluster_index)==0):
            print("ERROR: Cluster", GCname,"not found in", fname_baumgardt)
        elif (sum(cluster_index) > 1):
            print("ERROR: Cluster", GCname,"Found multiple times in", fname_baumgardt)
        else:
            half_mass_radius = fits_baumgardt[1].data['rh_m'][cluster_index][0]
            a = ((1/2)**(-2/3) -1)**(1/2) * half_mass_radius / 1000 # characteristic radius in kpc
            Mtot = fits_baumgardt[1].data['Mass'][cluster_index][0] * 10**-7 # total mass                                
        MGCparam=[Mtot,a]

        allGCparam.append(MGCparam)

    fits_baumgardt.close()
    return allGCparam

def createfile(astrometry_file,gcparameters_file,filename):

    file1 = fits.open(astrometry_file)
    cluster_index=file1[1].data['Name'].replace(" ", "")

    file2 = fits.open(gcparameters_file)
    cluster_index2=file2[1].data['ID'].replace(" ", "")

    GClist = list(set(cluster_index) & set(cluster_index2))

    with open(filename, 'w') as f:
        for item in GClist:
            f.write("%s\n" % item)

def getALL(filename,gc_frame):
    
    allGCname=[]
    f = open(filename, 'r')
    for line in f :
        GCname = line.strip()
        allGCname.append(GCname)
    

    allxGC,allyGC,allzGC,allvxGC,allvyGC,allvzGC=allGCcoord(allGCname,gc_frame,0)
    allGCparam = allPlummer(allGCname)

    return allGCname,allxGC,allyGC,allzGC,allvxGC,allvyGC,allvzGC,allGCparam
    
def allGCcoord(allGCname,gc_frame,ierror):

    allxGC=[]
    allyGC=[]
    allzGC=[]
    allvxGC=[]
    allvyGC=[]
    allvzGC=[]
    
    """ ASSIGN RA,DEC,PARALLAX,RADIAL VELOCITY AND PROPER MOTIONS OF THE GLOBULAR CLUSTER """
    fits_vasiliev = fits.open(fname_vasiliev)
    for i in range(0,len(allGCname)):
        GCname = allGCname[i]
        cluster_index=fits_vasiliev[1].data['Name'].replace(" ", "")==GCname
        if(sum(cluster_index)==1):

            RA = fits_vasiliev[1].data['RAJ2000'][cluster_index][0] # degrees
            DEC= fits_vasiliev[1].data['DEJ2000'][cluster_index][0] # degrees
            D = fits_vasiliev[1].data['Dist'][cluster_index][0] # kpc
            PMRA = fits_vasiliev[1].data['pmRA'][cluster_index][0] # mas/year
            PMDEC = fits_vasiliev[1].data['pmDE'][cluster_index][0] # mas/year
            RV =  fits_vasiliev[1].data['HRV'][cluster_index][0] # km/s
            RV_err =  fits_vasiliev[1].data['e_HRV'][cluster_index][0] # km/s
            PMRA_err = fits_vasiliev[1].data['e_pmRA'][cluster_index][0] # mas/year
            PMDEC_err = fits_vasiliev[1].data['e_pmDE'][cluster_index][0] # mas/year
            D_err=0.046*D # see Vasiliev 2019

            if(ierror==1):
                D1 = np.random.normal(D, D_err, 1)
                PMRA1 = np.random.normal(PMRA, PMRA_err, 1)
                PMDEC1 = np.random.normal(PMDEC, PMDEC_err, 1)
                RV1 = np.random.normal(RV, RV_err, 1)
                D=D1
                PMRA=PMRA1
                PMDEC=PMDEC1
                RV=RV1
            #print(D,PMRA,PMDEC,RV)
            #c1 = coord.SkyCoord(ra=RA*u.degree, dec=DEC*u.degree,distance=(PLX*u.mas).to(u.pc, u.parallax()),pm_ra_cosdec=PMRA*u.mas/u.yr,pm_dec=PMDEC*u.mas/u.yr,radial_velocity=RV*u.km/u.s,frame='icrs')
            c1 = coord.SkyCoord(ra=RA*u.degree, dec=DEC*u.degree,distance=D*1000*u.pc,pm_ra_cosdec=PMRA*u.mas/u.yr,pm_dec=PMDEC*u.mas/u.yr,radial_velocity=RV*u.km/u.s,frame='icrs')
        else:
            print("Warning, coordinates not found for", GCname, "in the file:",fname_vasiliev)

        """ AND TRANSFORM THEM TO GALACTOCENTRIC COORDINATES """
        gc2 = c1.transform_to(gc_frame)
    
        xGC=gc2.x.value/1000. # in units of kpc                                                                                               
        yGC=gc2.y.value/1000.
        zGC=gc2.z.value/1000.
        vxGC=gc2.v_x.value/10.# in units of 10km/s                                                                                            
        vyGC=gc2.v_y.value/10.
        vzGC=gc2.v_z.value/10.

        allxGC.append(xGC)
        allyGC.append(yGC)
        allzGC.append(zGC)
        allvxGC.append(vxGC)
        allvyGC.append(vyGC)
        allvzGC.append(vzGC)
        
    fits_vasiliev.close()
    return allxGC,allyGC,allzGC,allvxGC,allvyGC,allvzGC


def GCcoord(GCname,gc_frame,ierror):
    """ ASSIGN RA,DEC,PARALLAX,RADIAL VELOCITY AND PROPER MOTIONS OF THE GLOBULAR CLUSTER """
    fits_vasiliev = fits.open(fname_vasiliev)
    cluster_index=fits_vasiliev[1].data['Name'].replace(" ", "")==GCname
    if(sum(cluster_index)==1):

        RA = fits_vasiliev[1].data['RAJ2000'][cluster_index][0] # degrees                                                                                  
        DEC= fits_vasiliev[1].data['DEJ2000'][cluster_index][0] # degrees                                                                                  
        D = fits_vasiliev[1].data['Dist'][cluster_index][0] # kpc                                                                                          
        PMRA = fits_vasiliev[1].data['pmRA'][cluster_index][0] # mas/year                                                                                  
        PMDEC = fits_vasiliev[1].data['pmDE'][cluster_index][0] # mas/year                                                                                 
        RV =  fits_vasiliev[1].data['HRV'][cluster_index][0] # km/s                                                                                        
        RV_err =  fits_vasiliev[1].data['e_HRV'][cluster_index][0] # km/s                                                                                  
        PMRA_err = fits_vasiliev[1].data['e_pmRA'][cluster_index][0] # mas/year                                                                            
        PMDEC_err = fits_vasiliev[1].data['e_pmDE'][cluster_index][0] # mas/year                                                                           
        D_err=0.046*D # see Vasiliev 2019                                                                                                                  

        if(ierror==1):
            D1 = np.random.normal(D, D_err, 1)
            PMRA1 = np.random.normal(PMRA, PMRA_err, 1)
            PMDEC1 = np.random.normal(PMDEC, PMDEC_err, 1)
            RV1 = np.random.normal(RV, RV_err, 1)
            D=D1
            PMRA=PMRA1
            PMDEC=PMDEC1
            RV=RV1
        print(D,PMRA,PMDEC,RV)
        #c1 = coord.SkyCoord(ra=RA*u.degree, dec=DEC*u.degree,distance=(PLX*u.mas).to(u.pc, u.parallax()),pm_ra_cosdec=PMRA*u.mas/u.yr,pm_dec=PMDEC*u.mas/u.yr,radial_velocity=RV*u.km/u.s,frame='icrs')                                                                                                               
        c1 = coord.SkyCoord(ra=RA*u.degree, dec=DEC*u.degree,distance=D*1000*u.pc,pm_ra_cosdec=PMRA*u.mas/u.yr,pm_dec=PMDEC*u.mas/u.yr,radial_velocity=RV*u.km/u.s,frame='icrs')
    else:
        print("Warning, coordinates not found for", GCname, "in the file:",fname_vasiliev)

    """ AND TRANSFORM THEM TO GALACTOCENTRIC COORDINATES """
    gc2 = c1.transform_to(gc_frame)

    xGC=gc2.x.value/1000. # in units of kpc                                                                                                                
    yGC=gc2.y.value/1000.
    zGC=gc2.z.value/1000.
    vxGC=gc2.v_x.value/10.# in units of 10km/s                                                                                                             
    vyGC=gc2.v_y.value/10.
    vzGC=gc2.v_z.value/10.

    fits_vasiliev.close()
    return xGC,yGC,zGC,vxGC,vyGC,vzGC


def Nbody(GCname):

    """                                                                                                                                  READ INITIAL CONDITIONS                                                                                                              """

    filename="../plummer_profiles/"+GCname+'Plummer.fits'
    orbit=fits.open(filename)
    tbdata=orbit[1].data
    x0=tbdata.field('x0')
    y0=tbdata.field('y0')
    z0=tbdata.field('z0')
    vx0=tbdata.field('vx0')
    vy0=tbdata.field('vy0')
    vz0=tbdata.field('vz0')

    return x0,y0,z0,vx0,vy0,vz0



import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import astropy.coordinates as coord
import astropy.units as u
import inputMW



namefile='output126.bin'
f = FortranFile(namefile, 'r')
x = f.read_reals(dtype='float32')
y = f.read_reals(dtype='float32')
z = f.read_reals(dtype='float32')
vx = f.read_reals(dtype='float32')
vy = f.read_reals(dtype='float32')
vz = f.read_reals(dtype='float32')
f.close() 

print(len(x))

plt.figure()
plt.scatter(x,y)
plt.show()

gc_frame=inputMW.ref_frame()

c1 = coord.SkyCoord(x=x * u.kpc, y=y * u.kpc,z=z * u.kpc, v_x=vx*10 * u.km/u.s, v_y=vy*10 * u.km/u.s, v_z=vz*10 * u.km/u.s, frame=gc_frame)
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


plt.figure()
plt.scatter(l,b)
plt.show()

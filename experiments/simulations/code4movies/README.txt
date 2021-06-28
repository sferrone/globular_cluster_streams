This creates binary files saving the position of each star in the globular cluster

(1) creates the files integrate.so. This way python can pass arguments to our fourtran integrator. We do this in (2). When you choose your globular cluster, you have to do again edit the file GClistered (3). We need to have the plummer density profile created to run tpbody. If you want to do a new plummer profile, you need to do (4). You need to edit the code to point to the appropriate text file containing the names of the globular clusters. 

(1) f2py3 -c orbits.f90 -m integrate
(2) python3 TPBody.py
(3) edit the file:GClistered with names consistent with baumgard.fits
(4) python3 plummerIC.py
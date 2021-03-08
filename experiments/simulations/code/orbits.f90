MODULE integrate

REAL*8,DIMENSION(:),ALLOCATABLE :: xp,yp,zp,vxp,vyp,vzp,tesc,phiMW,phiGC

contains

  SUBROUTINE orbits(x,y,z,vx,vy,vz,GCparam,MWparam,model,modelshort,GCname,backward,pathorbits,deltat,Nstep,N)

!! Integrate the orbits of stars in three axisymmetric Galactic potentials:
!! described in Allen & Santillan 1991, Pouliasis et al 2017
!! =========================================================
!! UNITS USED 
!! Positions in units of kpc, velocities in units of 10km/s
!! time in units of 0.1 Gyr, as in Allen & Santillan 1991,Pouliasis et al 2017


  IMPLICIT NONE


  INTEGER,intent(in):: N,Nstep
  REAL,DIMENSION(N),INTENT(IN) :: x,y,z,vx,vy,vz
  REAL,DIMENSION(2),INTENT(IN) :: GCparam
  REAL,DIMENSION(18),INTENT(IN) :: MWparam
  REAL*8,INTENT(IN):: deltat
  CHARACTER*100,INTENT(IN) :: model,modelshort,backward,pathorbits,GCname

  REAL*8 :: deltat2
  
  INTEGER :: i,im,imm, istep
  INTEGER,DIMENSION(N) :: iesc
  REAL*8,DIMENSION(N) :: xnew,ynew,znew,vxnew,vynew,vznew,axnew,aynew,aznew,axp,ayp,azp
  REAL*8,DIMENSION(Nstep+1) :: xGCv,yGCv,zGCv,vxGCv,vyGCv,vzGCv
  REAL*8 :: xGC,yGC,zGC,vxGC,vyGC,vzGC,ekin,etot

  REAL*8 :: time,theta,theta0,omega_bar,energy
  CHARACTER*60 :: ll

  print*,'N,Nstep,deltat'
  print*,N,Nstep,deltat
  
 

  CALL ALLOCATE(N)

  CALL INPUTDATA(backward,model,modelshort,pathorbits,GCname,Nstep,omega_bar,xGCv,yGCv,zGCv,&
        vxGCv,vyGCv,vzGCv)

  IF(backward.EQ.'YES')&
       OPEN(UNIT=50,FILE='../outputDATA/'//TRIM(model)//TRIM(pathorbits)//'/orbit'//TRIM(GCname)//'.dat',FORM='FORMATTED')

  xGC=xGCv(Nstep+1)
  yGC=yGCv(Nstep+1)
  zGC=zGCv(Nstep+1)


  xp=x
  yp=y
  zp=z
  vxp=-vx
  vyp=-vy
  vzp=-vz
  
  theta0=-25.*atan(1.)/45. ! the bar is oriented of -25 degrees with respect to the x-y axis                                                                       
  IF(backward.EQ.'NO') THEN ! when the GC is treated as a N-body system, the initial bar angle is that found by rotating the current bar angle (theta0) by omega_bar* Nstep*deltat and taking into account that the rotation has imm = -1 
     imm = -1
     theta=theta0-imm*omega_bar*Nstep*deltat
     theta0=theta
  ENDIF

  im=1
  IF(backward.EQ.'YES') im=-1
  

  DO i=1,N
     CALL ACCEL(backward,xp(i),yp(i),zp(i),theta0,GCparam,MWparam,xGC,yGC,zGC,axp(i),ayp(i),azp(i),&
          phiMW(i),phiGC(i))
  ENDDO



  IF(backward.EQ.'YES')&
       WRITE(50,*)time,SNGL(xp),SNGL(yp),SNGL(zp),SNGL(vxp),SNGL(vyp),SNGL(vzp)
     

  deltat2=deltat*deltat

  tesc = -9999

  DO istep=1,Nstep
     time=istep*deltat
     theta=theta0-im*omega_bar*time ! attention, the sign of im needs to be checked yet


     xGC=xGCv(Nstep+1-istep)
     yGC=yGCv(Nstep+1-istep)
     zGC=zGCv(Nstep+1-istep)
     vxGC=-vxGCv(Nstep+1-istep)
     vyGC=-vyGCv(Nstep+1-istep)
     vzGC=-vzGCv(Nstep+1-istep)


     ! advance pos                                                                                                            

     xnew=xp+vxp*deltat+0.5*axp*deltat2
     ynew=yp+vyp*deltat+0.5*ayp*deltat2
     znew=zp+vzp*deltat+0.5*azp*deltat2
     ! calculate accelerations on new positions                                                                                

     DO i=1,N
        CALL ACCEL(backward,xnew(i),ynew(i),znew(i),theta,GCparam,MWparam,xGC,yGC,zGC,axnew(i),&
             aynew(i),aznew(i),phiMW(i),phiGC(i))
     ENDDO

     ! advance vel                                                                                                             
     vxnew=vxp+(axp+axnew)*deltat*0.5
     vynew=vyp+(ayp+aynew)*deltat*0.5
     vznew=vzp+(azp+aznew)*deltat*0.5
           
     xp=xnew
     yp=ynew
     zp=znew
     vxp=vxnew
     vyp=vynew
     vzp=vznew
     axp=axnew
     ayp=aynew
     azp=aznew

          ! estimate energy relative to the GC position and motion
     IF(backward.EQ.'NO') THEN
      DO I=1,N
         ekin=0.5*((vxp(i)-vxGC)**2+(vyp(i)-vyGC)**2+(vzp(i)-vzGC)**2)
         etot = ekin+phiGC(i)-phiMW(i)
         IF ((etot > 0).AND.(iesc(i).EQ.0)) THEN
            tesc(i)=time
            iesc(i)=1
         ENDIF
      ENDDO
   ENDIF


     IF(backward.EQ.'YES')&
          WRITE(50,*)time,SNGL(xnew),SNGL(ynew),SNGL(znew),SNGL(vxnew),SNGL(vynew),SNGL(vznew)


  ENDDO


  IF(backward.EQ.'YES')&
       CLOSE(UNIT=50)


END SUBROUTINE orbits



SUBROUTINE ACCEL(backward,xsinp,ysinp,zsinp,thetainp,GCparam,MWparam,xGC,yGC,zGC,axsout,aysout,azsout,phis,phiGC)

  IMPLICIT NONE

  REAL*8,INTENT(IN) :: xsinp,ysinp,zsinp,thetainp
  REAL,DIMENSION(2),INTENT(IN):: GCparam
  REAL,DIMENSION(18),INTENT(IN) :: MWparam
  CHARACTER*100,INTENT(IN) :: backward

  REAL*8,INTENT(OUT) :: axsout,aysout,azsout,phis,phiGC


  REAL*8 :: xs,ys,zs,xs2,ys2,zs2,rs,rs2,rs3,zmod,zmod2,xp,yp,zp
  REAL*8 :: Mbulge_r,Mhalo_r
  REAL*8 :: Md_thinmod,Md_intermod,Md_thickmod
  REAL*8 :: factbulge,facthalo,factthin,factthick,fact1,fact2,Tplus,Tminus,abarx, abary, abarz,ahalox,ahaloy,ahaloz
  REAL*8,DIMENSION(3) :: abulge,ahalo,athin,athick,abar,accGC
  REAL*8 :: phibulge,phihalo,phithin,phithick,phibar
  REAL*8 :: Md_bulge,rd_bulge,zd_bulge
  REAL*8 :: Md_halo,rd_halo,zd_halo,axhalo,bxhalo,cxhalo,thetahalo
  REAL*8 :: Md_thin,rd_thin,zd_thin
  REAL*8 :: Md_thick,rd_thick,zd_thick
  REAL*8 :: Md_bar,axbar,bxbar,cxbar
  REAL*8 :: MGC,aGC,xGC,yGC,zGC,rgc,rgc2,rgc3,Mr,amod

  REAL,PARAMETER :: G=1.

  xs=xsinp
  ys=ysinp
  zs=zsinp
  xs2=xs*xs
  ys2=ys*ys
  zs2=zs*zs
  rs=sqrt(xs2+ys2+zs2)
  rs2=rs*rs
  rs3=rs2*rs

  abulge=0.
  athin=0.
  athick=0.
  ahalo=0.

  MGC=GCparam(1)
  aGC=GCparam(2)


  Md_bulge=MWparam(1)
  Md_halo=MWparam(2)
  Md_thin=MWparam(3)
  Md_thick=MWparam(4)
  Md_bar=MWparam(5)
  zd_bulge=MWparam(6)
  rd_bulge=MWparam(7)
  zd_halo=MWparam(8)
  axhalo=MWparam(8)
  bxhalo=MWparam(9)
  cxhalo=MWparam(10)
  thetahalo=MWparam(11)
  zd_thin=MWparam(12)
  rd_thin=MWparam(13)
  zd_thick=MWparam(14)
  rd_thick=MWparam(15)
  axbar=MWparam(16)
  bxbar=MWparam(17)
  cxbar=MWparam(18)

  ! Miyamoto-Nagai/Plummer bulge                                                                                                        
  zmod=rd_bulge+(zs2+zd_bulge*zd_bulge)**0.5
  zmod2=zmod*zmod
  factbulge=-g*Md_bulge/(xs2+ys2+zmod2)**1.5
  abulge(1)=factbulge*xs
  abulge(2)=factbulge*ys
  abulge(3)=factbulge*zs*zmod/sqrt(zs2+zd_bulge*zd_bulge)
  phibulge=-g*Md_bulge/(xs2+ys2+zmod2)**0.5


  IF((bxhalo.EQ.0).AND.(cxhalo.EQ.0)) THEN
    zd_halo = axhalo
    ! Allen and Santillan halo
    Mhalo_r=Md_halo*(rs/zd_halo)**2.02
    Mhalo_r=Mhalo_r/(1+(rs/zd_halo)**1.02)
    facthalo=-g*Mhalo_r/rs3
    ahalo(1)=facthalo*xs
    ahalo(2)=facthalo*ys
    ahalo(3)=facthalo*zs
    fact1=-1.02/(1+(100./zd_halo)**1.02)+log(1+(100./zd_halo)**1.02)
    fact2=-1.02/(1+(rs/zd_halo)**1.02)+log(1+(rs/zd_halo)**1.02)
    phihalo=-(g*Mhalo_r/rs)-g*Md_halo/1.02/zd_halo*(fact1-fact2)
 ENDIF

 IF((bxhalo.GT.0).OR.(cxhalo.GT.0)) THEN
  ! Triaxial halo
  xp=xs*cos(thetahalo)+ys*sin(thetahalo)
  yp=-xs*sin(thetahalo)+ys*cos(thetahalo)
  zp=zs
  Tplus=sqrt((axhalo+xp)**2.+yp*yp+(bxhalo+sqrt(cxhalo*cxhalo+zp*zp))**2.)
  Tminus=sqrt((axhalo-xp)**2.+yp*yp+(bxhalo+sqrt(cxhalo*cxhalo+zp*zp))**2.)
  phihalo=-(g*Md_halo/2./axhalo)*log((xp-axhalo+Tminus)/(xp+axhalo+Tplus))
  ahalox=-2.*g*Md_bar*xp/((Tplus*Tminus)*(Tplus+Tminus))
  ahaloy=-g*Md_halo*yp/((2.*Tplus*Tminus)*(yp*yp+(bxhalo+sqrt(zp*zp+cxhalo*cxhalo))**2.))*(Tplus+Tminus-4*xp*xp/(Tplus+Tminus))
  ahaloz=-g*Md_halo*zp/((2.*Tplus*Tminus)*(yp*yp+(bxhalo+sqrt(zp*zp+cxhalo*cxhalo))**2.))*(Tplus+Tminus-4*xp*xp/(Tplus+Tminus))*&
    ((bxhalo+sqrt(zp*zp+cxhalo*cxhalo))/sqrt(zp*zp+cxhalo*cxhalo))
  ahalo(1)=ahalox*cos(-thetahalo)+ahaloy*sin(-thetahalo)
  ahalo(2)=-ahalox*sin(-thetahalo)+ahaloy*cos(-thetahalo)
  ahalo(3)=ahaloz
ENDIF

  ! Miyamoto-Nagai thin disk                                                                                                            
  Md_thinmod=Md_thin
  zmod=rd_thin+(zs2+zd_thin*zd_thin)**0.5
  zmod2=zmod*zmod
  factthin=-g*Md_thinmod/(xs2+ys2+zmod2)**1.5
  athin(1)=factthin*xs
  athin(2)=factthin*ys
  athin(3)=factthin*zs*zmod/sqrt(zs2+zd_thin*zd_thin)
  phithin=-g*Md_thinmod/(xs2+ys2+zmod2)**0.5

  ! Miyamoto-Nagai thick disk                                                                                                           
  Md_thickmod=Md_thick
  zmod=rd_thick+(zs2+zd_thick*zd_thick)**0.5
  zmod2=zmod*zmod
  factthick=-g*Md_thickmod/(xs2+ys2+zmod2)**1.5
  athick(1)=factthick*xs
  athick(2)=factthick*ys
  athick(3)=factthick*zs*zmod/sqrt(zs2+zd_thick*zd_thick)
  phithick=-g*Md_thickmod/(xs2+ys2+zmod2)**0.5

  abar(1)=0.
  abar(2)=0.
  abar(3)=0.
  phibar=0.
  

  IF(Md_bar.GT.0) THEN
     ! Stellar bar, Long & Murali 1992
     xp=xs*cos(thetainp)+ys*sin(thetainp)
     yp=-xs*sin(thetainp)+ys*cos(thetainp)
     zp=zs
     Tplus=sqrt((axbar+xp)**2.+yp*yp+(bxbar+sqrt(cxbar*cxbar+zp*zp))**2.)
     Tminus=sqrt((axbar-xp)**2.+yp*yp+(bxbar+sqrt(cxbar*cxbar+zp*zp))**2.)
     phibar=-(g*Md_bar/2./axbar)*log((xp-axbar+Tminus)/(xp+axbar+Tplus))
     abarx=-2.*g*Md_bar*xp/((Tplus*Tminus)*(Tplus+Tminus))
     abary=-g*Md_bar*yp/((2.*Tplus*Tminus)*(yp*yp+(bxbar+sqrt(zp*zp+cxbar*cxbar))**2.))*(Tplus+Tminus-4*xp*xp/(Tplus+Tminus))
     abarz=-g*Md_bar*zp/((2.*Tplus*Tminus)*(yp*yp+(bxbar+sqrt(zp*zp+cxbar*cxbar))**2.))*(Tplus+Tminus-4*xp*xp/(Tplus+Tminus))*&
       ((bxbar+sqrt(zp*zp+cxbar*cxbar))/sqrt(zp*zp+cxbar*cxbar))
     abar(1)=abarx*cos(-thetainp)+abary*sin(-thetainp)
     abar(2)=-abarx*sin(-thetainp)+abary*cos(-thetainp)
     abar(3)=abarz
  ENDIF

  accGC=0.
  phiGC=0.

  IF(backward.EQ.'NO') THEN ! add GC potential
     rgc=sqrt((xs-xGC)**2+(ys-yGC)**2+(zs-zGC)**2)
     rgc2=rgc*rgc
     rgc3=rgc2*rgc
     Mr=MGC*rgc*rgc*rgc/(rgc**2+aGC**2)**1.5 ! mass inside distance r from the centre of the Plummer sphere
     amod=-g*Mr/rgc3
     accGC(1)=amod*(xs-xGC)
     accGC(2)=amod*(ys-yGC)
     accGC(3)=amod*(zs-zGC)
     phiGC=-g*MGC/sqrt(rgc**2+aGC**2)
  ENDIF


  axsout=abulge(1)+ahalo(1)+athin(1)+athick(1)+abar(1)+accGC(1)
  aysout=abulge(2)+ahalo(2)+athin(2)+athick(2)+abar(2)+accGC(2)
  azsout=abulge(3)+ahalo(3)+athin(3)+athick(3)+abar(3)+accGC(3)

  phis=phibulge+phihalo+phithin+phithick+phibar ! this is the gravitational potential due to the Galaxy only. The GC potential is saved apart, as phiGC


END SUBROUTINE ACCEL


SUBROUTINE ALLOCATE(N)

  IMPLICIT NONE

  INTEGER :: N

  ALLOCATE(xp(N))
  ALLOCATE(yp(N))
  ALLOCATE(zp(N))
  ALLOCATE(vxp(N))
  ALLOCATE(vyp(N))
  ALLOCATE(vzp(N))
!  ALLOCATE(axp(N))
!  ALLOCATE(ayp(N))
!  ALLOCATE(azp(N))
  ALLOCATE(tesc(N))
  ALLOCATE(phiMW(N))
  ALLOCATE(phiGC(N))


END SUBROUTINE ALLOCATE

SUBROUTINE INPUTDATA(backward,model,modelshort,pathorbits,GCname,Nstep,omega_bar,xGCv,yGCv,zGCv,&
  vxGCv,vyGCv,vzGCv)

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: Nstep
  CHARACTER*100,INTENT(IN) :: model,modelshort,backward,GCname, pathorbits
  REAL*8,INTENT(OUT) :: omega_bar
  REAL*8,DIMENSION(Nstep+1),INTENT(OUT) :: xGCv,yGCv,zGCv,vxGCv,vyGCv,vzGCv


  INTEGER :: istep
  REAL*8 :: time


  OPEN(UNIT=30,FILE='../inputDATA/'//TRIM(modelshort)//'_barpatternspeed')
  READ(30,*)omega_bar
  print*,'omegabar: ',omega_bar
  CLOSE(UNIT=30)

  xGCv=0.
  yGCv=0.
  zGCv=0.
  vxGCv=0.
  vyGCv=0.
  vzGCv=0.

  IF(backward.EQ.'NO') THEN
     OPEN(UNIT=50,FILE='../outputDATA/'//TRIM(model)//TRIM(pathorbits)//'/orbit'//TRIM(GCname)//'.dat',FORM='FORMATTED')
     DO istep=1,Nstep+1
      READ(50,*)time,xGCv(istep),yGCv(istep),zGCv(istep),vxGCv(istep),vyGCv(istep),vzGCv(istep)
     ENDDO
  ENDIF


END SUBROUTINE INPUTDATA

SUBROUTINE DEALLOCATE

  IMPLICIT NONE


  DEALLOCATE(xp)
  DEALLOCATE(yp)
  DEALLOCATE(zp)
  DEALLOCATE(vxp)
  DEALLOCATE(vyp)
  DEALLOCATE(vzp)
!  DEALLOCATE(axp)
!  DEALLOCATE(ayp)
!  DEALLOCATE(azp)
  DEALLOCATE(tesc)
  DEALLOCATE(phiMW)
  DEALLOCATE(phiGC)

END SUBROUTINE DEALLOCATE


END MODULE INTEGRATE



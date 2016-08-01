c ********************
      program genmod
c
c  $Header$
c  this program finds the complex love nos for Europa, assuming a Maxwell
c     solid ice layer.
c
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      dimension r(numb),rho(numb),vp(numb),vs(numb),g(numb)
     2  ,vfsq(numb),kdis(10)
      common /cgen/ r,rho,vp,vs
     2  ,vfsq,ndis,kdis,npts,qrhobar,z
      real rhobar, period,visc,viscbottom,vpc,vsc,thickf,densf,vpf,
     2  radius,densi, thradc, thbottom, Erock, nurock
      real rhorock, murock, lambdarock, Eice, nuice, muice, lambdaice
c     
c  file 10 has input values
      open (10, file='in.love',status='unknown')
c  the model at the node points will be written on file 16
      open (16, file='out.model',status='unknown')
c  the input values and the Love numbers will be written on file 20
      open (20, file='out.love',status='unknown')
c
cccccccccccccccccccccccccccccccccccccccccccc
c
c   rhobar is the mean density of Europa (cgs)
c     rhobar=3.040
      read (10,*) rhobar
      write (20,*) rhobar
542   format(f20.10)
c
c  when maxwell=1, let the ice be a Maxwell solid
c  when maxwell=0, let the ice be elastic
c     maxwell=1
      read (10,*) maxwell
      write (20,*) maxwell
c  period is the tidal period, in days
c     period=3.5
      read (10,*) period
      write (20,*) period
c  convert to seconds
      period=period*86400.
      freq=2.*3.14159/period
      ci=(0.,1.)
c  visc is the viscosity of the ice, in Pa-sec
c     visc=1.e21
c  viscbottom is the viscosity of the bottom ice, Pa-sec
      read (10, *) visc
      write (20,544) visc
544   format(e12.5)
c  read in the bottom viscosity from the input file.
      read (10, *) viscbottom
      write (20,512) viscbottom
512   format(e12.5)
c  convert to cgs
      visc=visc*10.0
      viscbottom=viscbottom*10.0
c
c  radc is the outer radius of the core (in km)
c  vpc is the p-wave velocity in the core (in km/s)
c  vsc is the s-wave velocity in the core (in km/s)
c     vpc=11.0     
c     vsc=6.0
c      read (10,*) vpc
c      write (20,*) vpc
c      read (10,*) vsc
c      write (20,*) vsc
c Changed so we input Young's Mod, poisson's ratio and
c rock density rather than p and s wave velocities.
c Rock density is in cgs units.  E is in Pa
      read (10, *) Erock
      write (20, *) Erock
      read (10, *) nurock
      write (20, *) nurock
      read (10, *) rhorock
      write (20, *) rhorock
c convert densities from cgs to mks units      
      rhorock=1000.0*rhorock
      
      murock=Erock/(2.0*(1.0+nurock))
      lambdarock=Erock*nurock/((1.0+nurock)*(1.0-(2.0*nurock)))
c program wants seismic velocities in cm/s
      vpc=100.0*sqrt((2.0*murock+lambdarock)/(rhorock))
      vsc=100.0*sqrt(murock/rhorock)
c      write(20,544) vpc
c      write(20,544) vsc      
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Amy and McCall changed input format 4.10.2005      ccc
c      cccccccccccccccccccccccccccccccccccccccccccccccccc
c Changed so we input Young's Mod, poisson's ratio and
c rock density rather than p and s wave velocities.
c Rock density is in cgs units.  E is in Pa
      read (10, *) Eice
      write (20, *) Eice
      read (10, *) nuice
      write (20, *) nuice
      read (10, *) densi
      write (20, *) densi
c convert densities from cgs to mks units      
      densi=1000.0*densi
      
      muice=Eice/(2.0*(1.0+nuice))
      lambdaice=Eice*nuice/((1.0+nuice)*(1.0-(2.0*nuice)))
c program wants seismic velocities in cm/s
      vpi=100.0*sqrt((2.0*muice+lambdaice)/(densi))
      vsi=100.0*sqrt(muice/densi)
c      write(20,544) muice
c      write(20,544) lambdaice
      
      densi=densi/1000.0

c  ifluid=1 if we want a liquid water layer beneath the ice.
c    ifluid=0 if we do not
c  thickf is the thickness of the fluid layer (in km)
c  densf is the density of the fluid layer (in cgs)
c  vpf is the p-wave velocity in the fluid layer (in km/s)
c     ifluid=1   
c     thickf=75.
c     densf=1.0     
c     vpf=1.45 
      read (10,*) ifluid
      write (20,*) ifluid
      read (10,*) thickf
      write (20,*) thickf
      read (10,*) densf
      write (20,*) densf
      read (10,*) vpf
      write (20,*) vpf
      thickf=thickf*1.e5
      vpf=vpf*1.e5
c
c  radius is the outer radius of the planet (in km)
c  densi is the density of the outer ice layer (in cgs)
c  vpii is the p-wave velocity in the outer ice layer (in km/s)
c  vsii is the s-wave velocity in the outer ice layer (in km/s)
c     radius=1569.   
c     densi=0.94     
c     vpii=2.5     
c     vsii=1.25
      read (10,*) radius
      write (20,*) radius
c      read (10,*) vpii
c      write (20,*) vpii
c      read (10,*) vsii
c      write (20,*) vsii
      radius=radius*1.e5
c      vpi=vpii*1.e5
c      vsi=vsii*1.e5
c
c  make seismic velocities complex
c  vsibottom, vpibottom store seismic velocities for bottom ice layer.
      if (maxwell.eq.1) then
      vsitemp=vsi
      vpitemp=vpi
      tau=visc/vsitemp**2/densi
      taubottom=viscbottom/vsitemp**2/densi
      vsi=vsitemp*sqrt(ci*freq/(ci*freq+1./tau))
      vpi=vpitemp*sqrt(1./(ci*freq+1./tau))*
     2  sqrt(ci*freq+1./tau*(1.-4./3.*(vsitemp/vpitemp)**2))

c  Calculate seismic velocities for the bottom layer.
      vsibottom=vsitemp*sqrt(ci*freq/(ci*freq+1./taubottom))
      vpibottom=vpitemp*sqrt(1./(ci*freq+1./taubottom))*
     2  sqrt(ci*freq+1./taubottom*(1.-4./3.*(vsitemp/vpitemp)**2))
      endif
c
c   if ifluid=0, then rfl1 and rfl2 are the radii of the solid
c     inner core and fluid outer core (in km).  Otherwise, these
c     are not used
      rfl1=10.0        
      rfl2=20.0
      rfl1=rfl1*1.e5
      rfl2=rfl2*1.e5
c
      if (real(rfl1).gt.real(rfl2)) then
      print *,'error #3',rfl1,rfl2
      stop
      endif
c
100   format(i5,3f10.0)
101   format(4f10.0)
c
      coef=4./3.*3.14159*6.67e-8
c
cccccccccccccccccccccccccccccccccccccc
c  added in one more input parameter to set the location of
c  the discontinuity within the ice.  AB
c  thradc=thickness of ice layer (in km)
c  thbottom=thickness of the bottom (presumably warmer)ice layer.  
c     thradc=2.0
      read (10,*) thradc
      write (20,*) thradc
      read (10,*) thbottom
      write (20, *) thbottom
      thradc=thradc
      thbottom=thbottom
c  thickif is the thickness of the combined ocean+ice layer
      thickif=thickf/1.e5+thradc+thbottom
c  radc is radius (in km) of rocky core.  It is not read in.
c    Instead, it is determined from all the other radii
      radc=radius/1.e5-thickif
      radc=radc*1.e5
      if (real(radc).gt.real(radius)) then
      print *,'error #1',radc,radius
      stop
      endif
      if (real(rfl2).ge.real(radc)) then
      print *,'error #4',rfl2,radc
      stop
      endif
      if ((ifluid.eq.1).and.((real(radc)+thickf).gt.radius)) then
      print *,'error #2',radc,thickf,radius
      stop
      endif
c
c  IF THERE'S AN OCEAN, DO THIS STUFF...
      if (ifluid.eq.1) then

      radf=radc+thickf
      densc=(rhobar*radius**3-densf*(radf**3-radc**3)-
     2     densi*(radius**3-radf**3))/radc**3

c  npts is the total no of nodes in the structural model
c  ndis is the no of discontinuities (should be 2 for a single
c    ice layer over a single fluid layer over a single rocky layer)
c  kdis(1) is the node number of the discontinuity between
c   the inner-most and next-inner-most layer
c  kdis(2) is the node number for the next discontinuity
c  etc
c  nominal values I tried out 1st time: 
c     npts=165
c     ndis=3
c     kdis(1)=120
c     kdis(2)=142
c     kdis(3)=150
      read (10,*) npts
      write (20,*) npts
      read (10,*) ndis
      write (20,*) ndis
      do 642 nn=1,ndis
      read (10,*) kdis(nn)
      write (20,*) kdis(nn)
642   continue
543   format(i5)
c  set parameters in solid core
      dr=radc/float(kdis(1))
      do 10 i=1,kdis(1)
      r(i)=dr*float(i)
      rho(i)=densc
      vp(i)=vpc
      vs(i)=vsc
      g(i)=coef*r(i)*densc
      vfsq(i)=-(g(i)/vpc)**2
10    continue
      r(kdis(1))=radc
c  set parameters in liquid water layer
      dr=0.
      do 11 i=kdis(1)+1,kdis(2)
      r(i)=r(i-1)+dr
      dr=thickf/float(kdis(2)-kdis(1)-1)
      rho(i)=densf
      vp(i)=vpf
      vs(i)=0.
      g(i)=coef/r(i)**2*
     2 (densc*radc**3+densf*(r(i)**3-radc**3))
      vfsq(i)=-(g(i)/vpf)**2
11    continue
      r(kdis(2))=radc+thickf
c  set parameters in bottom warm ice layer  AB
c  bottom and top ice layers will need to have different seismic velocities in general?
      dr=0.
      do 12 i=kdis(2)+1,kdis(3)
      r(i)=r(i-1)+dr
      dr=(thbottom*1.e5)/float(kdis(3)-kdis(2)-1)
      rho(i)=densi
      vp(i)=vpibottom
      vs(i)=vsibottom
      g(i)=coef/r(i)**2*
     2 (densc*radc**3+
     2 densf*(radf**3-radc**3)+
     2 densi*(r(i)**3-radf**3))
c    2 densf*(r(i)**3-radf**3))
      vfsq(i)=-(g(i)/vpibottom)**2
c     vfsq(i)=-(g(i)/vpi)**2
12    continue
      r(kdis(3))=radc+thickf+thbottom*1.e5
c   set parameters for outermost ice layer.  
      dr=0.
      do 13 i=kdis(3)+1,npts
      r(i)=r(i-1)+dr
      dr=(thradc*1.e5)/float(npts-kdis(3)-1)
      rho(i)=densi
      vp(i)=vpi
      vs(i)=vsi
      g(i)=coef/r(i)**2*
     2 (densc*radc**3+
     2 densf*(radf**3-radc**3)+
     2 densi*(r(i)**3-radf**3))
c    2 densf*(r(i)**3-radf**3))
      vfsq(i)=-(g(i)/vpi)**2
13    continue
      r(npts)=radius
      rhobar=rhobar
      
      
c  IF THERE IS NOT AN OCEAN, DO THIS STUFF... No provisions
c  made in here yet for a 2nd ice layer.  
      else

      densc=(rhobar*radius**3-
     2     densi*(radius**3-radc**3))/radc**3
c     print *,'densc=',real(densc)

c     ndis=3
c     kdis(1)=6
c     kdis(2)=12
c     kdis(3)=138
c     npts=165
c  set parameters in dummy solid core
      dr=rfl1/float(kdis(1))
      do 20 i=1,kdis(1)
      r(i)=dr*float(i)
      rho(i)=densc
      vp(i)=vpc
      vs(i)=vsc
      g(i)=coef*r(i)*densc
      vfsq(i)=-(g(i)/vpc)**2
20    continue
      r(kdis(1))=rfl1
c  set parameters in dummy liquid outer core
      dr=0.
      do 21 i=kdis(1)+1,kdis(2)
      r(i)=r(i-1)+dr
      dr=(rfl2-rfl1)/float(kdis(2)-kdis(1)-1)
      rho(i)=densc
      vp(i)=vpc
      vs(i)=0.
      g(i)=coef*r(i)*densc
      vfsq(i)=-(g(i)/vpc)**2
21    continue
      r(kdis(2))=rfl2
c  set parameters in real solid core
      dr=0.
      do 22 i=kdis(2)+1,kdis(3)
      r(i)=r(i-1)+dr
      dr=(radc-rfl2)/float(kdis(3)-kdis(2)-1)
      rho(i)=densc
      vp(i)=vpc
      vs(i)=vsc
      g(i)=coef*r(i)*densc
      vfsq(i)=-(g(i)/vpc)**2
22    continue
      r(kdis(3))=radc
c  set parameters in outer ice layer
      dr=0.
      do 23 i=kdis(3)+1,npts
      r(i)=r(i-1)+dr
      dr=(radius-radc)/float(npts-kdis(3)-1)
      rho(i)=densi
      vp(i)=vpi
      vs(i)=vsi
      g(i)=coef/r(i)**2*
     2 (densc*radc**3+densi*(r(i)**3-radc**3))
      vfsq(i)=-(g(i)/vpi)**2
23    continue
      r(npts)=radius
      rhobar=rhobar

      endif
c
      if ((kdis(ndis)+1).ge.npts) then
      print *,'error #6',ndis,kdis(ndis),npts
      stop
      endif
      if (npts.gt.numb) then
      print *,'error #5',npts,numb
      stop
      endif
c
      do 677 i=1,npts
      r(i)=r(i)/100.
      rho(i)=rho(i)*1000.
      vp(i)=vp(i)/100.
      vs(i)=vs(i)/100.
677   continue
c
      qrhobar=rhobar*1000.
      z=0.
      call main
c
      stop
      end

      subroutine MAIN
C   PROGRAM TO COMPUTE LOAD LOVE NUMBERS
C   THIS IS THE ORIGINAL COPY OF DAHLEN'S PROGRAM
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      dimension qr(numb),qrho(numb),qvp(numb),qvs(numb)
     2  ,qvfsq(numb),ikdis(10)
      common /cgen/ qr,qrho,qvp,qvs
     2  ,qvfsq,indis,ikdis,inpts,qrhobar,qz
      COMMON/E/NDIS,KDIS
      COMMON/F/C,NIR,NER
      COMMON/P/CAPW,EA,EAH,ELL
      COMMON/C/STEP
      COMMON/J/NPT
      COMMON/R2/R2
      DIMENSION C(numb),FLAM(numb),FMU(numb),G(numb),
     2  R(numb),R2(numb),RHO(numb)
     1,VFSQ(numb),VP(numb),VS(numb),BULK(numb),
     2  ELL(numb),STEP(8),KDIS(28),MOD
     2EL(5),LRD(36)
      PI=3.1415926536
      CON=PI*6.67D-11
      CAPW=7.292116D-05
      EA=(1.)/(2.982582)
      EPS=1.D-08
      NPT=1
      T3=2./3.
      n=inpts
      do 677 i=1,n
      r(i)=qr(i)
      rho(i)=qrho(i)
      vp(i)=qvp(i)
      vs(i)=qvs(i)
      vfsq(i)=qvfsq(i)
677   continue
      CALL DG
      NIR=1
      NER=N
      RN=R(N)
      rhobar=qrhobar
      z=qz
      RABOHR=RHOBAR
      ZED=.5*Z
      DO 620 I=1,N
      R(I)=R(I)/RN
      R2(I)=R(I)*R(I)
  620 RHO(I)=RHO(I)/RHOBAR
      CALL PREINT
  970 DO 971 I=1,N
  971 C(I)=RHO(I)
      CALL GLINT(Q1)
      DO 972 I=1,N
  972 C(I)=RHO(I)*R(I)*R(I)
      CALL GLINT(Q2)
c  the following 4 statements adjust the density so
c    that it gives the correct C/(M*r**2).  I want to comment them out.
c      A1=26.25*(ZED-Q2)-6.25+18.75*Q1
c      A2=43.75*(Q2-ZED)+8.75-26.25*Q1
c      DO 973 I=1,N
c  973 RHO(I)=RHO(I)-A1-A2*R(I)*R(I)
      DO 974 I=1,N
      R(I)=R(I)*RN
  974 RHO(I)=RHO(I)*RHOBAR
      CALL DG
      RHOBAR=RABOHR
      ndis=indis
      do 877 kk=1,ndis
      kdis(kk)=ikdis(kk)
877   continue
      GN=CON*RHOBAR*RN
      V=SQRT (GN*RN)
      WN=V/RN
      CAPW=CAPW/WN
C   CON IS VISCOUS PARAMETER.  ABOVE, IT WAS USED IN SUBROUTINE DG
      CON=1.
      DO 7 I=1,N
      R(I)=R(I)/RN
      VP(I)=VP(I)/V
      VS(I)=VS(I)/V
      RHO(I)=RHO(I)/RHOBAR
      FMU(I)=RHO(I)*VS(I)*VS(I)
      FLAM(I)=RHO(I)*VP(I)*VP(I)-2.*FMU(I)
      BULK(I)=FLAM(I)+T3*FMU(I)
    7 G(I)=G(I)/GN
      CALL ELLIP
      DO 701 I=1,N
      R(I)=R(I)*RN
      RHO(I)=RHO(I)*RHOBAR
      VP(I)=VP(I)*V
      VS(I)=VS(I)*V
  701 G(I)=G(I)*GN
   30 I=1
      do 2 i=1,n
      NP=NP+1
      if (i.eq.1) WRITE(16,902) NP,(MODEL(J),J=1,5),N
  902 FORMAT(1H1,116X,I3/10X,10HMODEL NAME,10X,6HLEVELS//5X,5A4,I10///4X
     1,5HLEVEL,2X,6HRADIUS,5X,3HRHO,6X,9HVP (REAL),4X,9hVP (IMAG),
     2  2x,9hVS (REAL),4x,9HVS (IMAG),3X,7HGRAVITY,3X,
     211HELLIPTICITY/)
    3 WRITE(16,9021) I,real(R(I)),real(RHO(I)),real(VP(I)),aimag(vp(i))
     2  ,real(VS(I)),aimag(vs(i)),real(G(I)),real(ELL(I))
 9021 FORMAT(4X,I3,F11.0,1f9.0,2x,2(F11.2,1x,f9.2,1x),1x,2E13.4)
2     continue
c      WRITE(16,903) RHOBAR,Z
  903 FORMAT(//10X,8HRHOBAR= ,F10.4,10X,3HZ= ,F10.8)
      EAR=1./EA
      EAHR=1./EAH
c      WRITE(16,9031) EAR,EAHR
 9031 FORMAT(  10X,23HELLIPTICITY = ONE OVER ,F7.3,/10X,35HHYDROSTATIC E
     1LLIPTICITY = ONE OVER ,F7.3)
      DO 702 I=1,N
      R(I)=R(I)/RN
      RHO(I)=RHO(I)/RHOBAR
      VP(I)=VP(I)/V
      VS(I)=VS(I)/V
  702 G(I)=G(I)/GN
      CALL STEPS(EPS,VERTNO)
      W=0.
      WSQ=W*W
c      WRITE(16,706)
  706 FORMAT(1H1,//////)
      lltot=1
      lrd(1)=2
      DO 705 LL=1,LLTOT
      L=LRD(LL)
      FL=L
      FL1=FL+1.
      FL2=FL+FL1
      FL3=FL*FL1
      SFL3=SQRT(FL3)
      CALL LOVENO(HPR,FKPR,FLPR)
  703 FORMAT(2HL=,I3,4x,2HH=,2(F12.8,1x),4x,2HK=,2(F12.8,1x),4x,2HL=
     1,2(F12.8,1x))
      WRITE(20,703) L,HPR,FKPR,FLPR
  705 CONTINUE
      return
      END
      SUBROUTINE LOVENO(HPR,FKPR,FLPR)
C     FOR USE WITH AN EARTH MODEL WITHOUT A SURFICIAL OCEANIC LAYER.
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/E/NDIS,KDIS
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),PHIC(numb),PHIPRC(numb),AS(6,3),B(3,3),BINV(3,3),AR(3),
     2ASA(6),KDIS(28)
      NIC=KDIS(1)
      NICP=NIC+1
      NOC=KDIS(2)
      NOCP=NOC+1
      IN=2
      IF(L.GT.20) IN=15
      CALL SPSJFG(IN,AS)
      DO 1 K=1,3
      AS(2,K)=AS(2,K)/SFL3
      AS(5,K)=-AS(5,K)
    1 AS(6,K)=-AS(6,K)
      NS=IN+1
      Y=R(IN)
      DO 4 I=NS,NIC
      X=Y
      Y=R(I)
    4 CALL SMTR(X,Y,AS)
      F1=AS(3,1)-RHO(NICP)*AS(5,1)-RHO(NICP)*G(NIC)*AS(1,1)
      F2=AS(3,2)-RHO(NICP)*AS(5,2)-RHO(NICP)*G(NIC)*AS(1,2)
      F3=AS(3,3)-RHO(NICP)*AS(5,3)-RHO(NICP)*G(NIC)*AS(1,3)
      Q1=AS(4,1)
      Q2=AS(4,2)
      Q3=AS(4,3)
      BIGA=-(F3/F2-Q3/Q2)/(F1/F2-Q1/Q2)
      BIGB=-(F3/F1-Q3/Q1)/(F2/F1-Q2/Q1)
      F=AS(5,3)+BIGB*AS(5,2)+BIGA*AS(5,1)
      FPR=AS(6,3)-(FL1/R(NIC))*AS(5,3)-4.*RHO(NICP)*AS(1,3)+BIGB*(AS(6
     1,2)-(FL1/R(NIC))*AS(5,2)-4.*RHO(NICP)*AS(1,2))+BIGA*(AS(6,1)-(FL
     21/R(NIC))*AS(5,1)-4.*RHO(NICP)*AS(1,1))
      PHIC(NICP)=F
      PHIPRC(NICP)=FPR
      NS=NICP+1
      Y=R(NICP)
      DO 30 I=NS,NOC
      X=Y
      Y=R(I)
      CALL CMTR(X,Y,F,FPR)
      PHIC(I)=F
   30 PHIPRC(I)=FPR
      ROC=RHO(NOC)
      GC=G(NOC)
      FACC=PHIPRC(NOC)/PHIC(NOC)+FL1/R(NOC)
      DO 5 J=1,6
      DO 5 K=1,3
    5 AS(J,K)=0.
      AS(1,1)=1.
      AS(3,1)=ROC*GC
      AS(6,1)=4.*ROC
      AS(2,2)=1.
      AS(5,3)=1.
      AS(3,3)=ROC
      AS(6,3)=FACC
      NS=NOCP+1
      Y=R(NOCP)
      DO 6 I=NS,N
      X=Y
      Y=R(I)
    6 CALL SMTR(X,Y,AS)
      DO 7 K=1,3
      B(1,K)=AS(3,K)
      B(2,K)=AS(4,K)
    7 B(3,K)=AS(6,K)
      CALL MATADJ(B,BINV)
      DET=B(1,1)*BINV(1,1)+B(1,2)*BINV(2,1)+B(1,3)*BINV(3,1)
      DO 48 J=1,3
      DO 48 K=1,3
   48 BINV(J,K)=BINV(J,K)/DET
      DO 47 J=1,3
c  The following is used for load love nos
c  47 AR(J)=-G(N)*BINV(J,1)-4.*BINV(J,3)
c  The following is used for body tide love nos
   47 AR(J)=-4.*BINV(J,3)
      DO 41 J=1,6
   41 ASA(J)=0.
      DO 42 J=1,6
      DO 42 K=1,3
   42 ASA(J)=ASA(J)+AS(J,K)*AR(K)
      FKPR=-1.-(FL2/(4.*R(N)))*ASA(5)
      HPR=G(N)*(FL2/(4.*R(N)))*ASA(1)
      FLPR=G(N)*(FL2/(4.*R(N)))*ASA(2)
      RETURN
      END
      SUBROUTINE DG
C     UNALTERED DPFG VERSION
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),X(4)
      G(1)=0.
      F=0.
      DO 10 I=2,N
      J=I-1
      RJI=R(I)-R(J)
    1 IF(real(RHO(I))) 6,6,2
    2 IF(real(R(I))-real(r(j))) 5,5,3
    3 DEL=RJI/3.
      X(1)=R(J)
      X(2)=X(1)+DEL
      X(4)=R(I)
      X(3)=X(4)-DEL
      DO 4 K=1,4
      DEL=(X(K)-R(J))/RJI
      RO=RHO(J)+DEL*(RHO(I)-RHO(J))
    4 X(K)=RO*X(K)*X(K)
      F=F+.125*RJI*(X(1)+X(4)+3.*(X(2)+X(3)))
    5 G(I)=4.*CON*F/(R(I)*R(I))
      GO TO 10
    6 IF(I.GT.2) GO TO 8
    7 ROP=0.
      GP=4.*CON*RHO(1)/3.
      GO TO 9
    8 C=VP(J)*VP(J)-4.*VS(J)*VS(J)/3.
      ROP=-RHO(J)*(VFSQ(J)/G(J)+G(J)/C)
      GP=4.*CON*RHO(J)-2.*G(J)/R(J)
    9 RHO(I)=RHO(J)
      G(I)=G(J)
      SR=.5*ROP
      SG=.5*GP
      ER=ROP
      EG=GP
      RHO(I)=RHO(I)+RJI*SR
      G(I)=G(I)+RJI*SG
      RIJ=.5*(R(I)+R(J))
      V=.5*(VP(I)+VP(J))
      S=.5*(VS(I)+VS(J))
      C=V*V-4.*S*S/3.
      VF=.5*(VFSQ(I)+VFSQ(J))
      ROP=-RHO(I)*(VF/G(I)+G(I)/C)
      GP=4.*CON*RHO(I)-2.*G(I)/RIJ
      A=.292893218814
      SR=A*(ROP-ER)
      SG=A*(GP-EG)
      ER=ER+3.*SR-A*ROP
      EG=EG+3.*SG-A*GP
      RHO(I)=RHO(I)+RJI*SR
      G(I)=G(I)+RJI*SG
      ROP=-RHO(I)*(VF/G(I)+G(I)/C)
      GP=4.*CON*RHO(I)-2.*G(I)/RIJ
      A=1.707106781186
      SR=A*(ROP-ER)
      SG=A*(GP-EG)
      ER=ER+3.*SR-A*ROP
      EG=EG+3.*SG-A*GP
      RHO(I)=RHO(I)+RJI*SR
      G(I)=G(I)+RJI*SG
      C=VP(I)*VP(I)-4.*VS(I)*VS(I)/3.
      ROP=-RHO(I)*(VFSQ(I)/G(I)+G(I)/C)
      GP=4.*CON*RHO(I)-2.*G(I)/R(I)
      A=1./6.
      SR=A*(ROP-2.*ER)
      SG=A*(GP-2.*EG)
      ER=ER+3.*SR-.5*ROP
      EG=EG+3.*SG-.5*GP
      RHO(I)=RHO(I)+RJI*(SR-ER/3.)
      G(I)=G(I)+RJI*(SG-EG/3.)
      F=G(I)*R(I)*R(I)/(4.*CON)
   10 CONTINUE
      RHOBAR=3.*F/(R(N)**3)
      RETURN
      END
      SUBROUTINE ELLIP
C     FOR USE IN LOVE NUMBER PROGRAMS.
C     DOES NOT DEFINE DELTA M SUB ELL AND DOES NOT CALL DERIV.
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/F/C,NIR,NER
      COMMON/P/CAPW,EA,EAH,ELL
      COMMON/MP/FMU1,BULK1,RHO1,G1
      COMMON/R2/R2
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),C(numb),R2(numb),P(numb),T(numb),ETA(numb),ELL(numb)
      NIR=1
      DO 1 I=1,N
    1 C(I)=RHO(I)
      NER=0
      DO 2 I=1,N
      NER=NER+1
      CALL GLINT(Q)
    2 P(I)=Q
      DO 3 I=1,N
    3 C(I)=RHO(I)*R2(I)
      NER=0
      DO 4 I=1,N
      NER=NER+1
      CALL GLINT(Q)
    4 T(I)=Q
      DO 5 I=2,N
      ZED=(2.*T(I))/(3.*R2(I)*P(I))
    5 ETA(I)=(25./4.)*(1.-1.5*ZED)*(1.-1.5*ZED)-1.
      ETA(1)=0.
      FM=.75*CAPW*CAPW
      EAH=(5.*FM)/(2.*(ETA(N)+2.))
      DO 6 I=2,N
    6 C(I)=ETA(I)/R(I)
      C(1)=0.
      NER=N
      CALL GLINT(Q)
      NER=0
      DO 7 I=1,N
      NER=NER+1
      CALL GLINT(Q1)
      ARG=Q1-Q
    7 ELL(I)=EAH*EXP(ARG)
      RETURN
      END
      SUBROUTINE GLINT(Q)
C     ALTERED DPFG VERSION  WORKS FOR ARBITRARY UPPER AND LOWER LIMITS
      implicit complex(a-h,o-z)
      parameter(numb=1000)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,O,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/F/F,NIR,NER
      COMMON/INT/B(numb)
      COMMON/R2/R2(numb)
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),F(numb),A(numb),C(numb)
      Q=0.
      IF(NER.EQ.NIR) RETURN
      NIRP=NIR+1
      NERM=NER-1
      IF(NERM.LT.NIRP) GO TO 3
      DO 1 I=NIRP,NERM
    1 Q=F(I)*B(I)+Q
    3 Q=Q+F(NIR)*A(NIR)+F(NER)*C(NER)
      RETURN
      ENTRY PREINT
      QTR=.25
      C6=1./6.
      C12=1./12.
      Y=0.
      YSQ=0.
      B(1)=0.
      DO 2 I=2,N
      IM1=I-1
      X=Y
      Y=R(I)
      H=Y-X
      XSQ=YSQ
      YSQ=R2(I)
      XY=C6*X*Y
      AI=(QTR*XSQ+XY+C12*YSQ)*H
      A(I)=AI
      B(IM1)=AI+B(IM1)
      C(I)=(QTR*YSQ+XY+C12*XSQ)*H
    2 B(I)=C(I)
      A(1)=B(1)
      RETURN
      END
      SUBROUTINE SHANKS
C     UNALTERED DPFG VERSION
      implicit complex(a-h,o-z)
      COMMON/I/B,C,DX,I
      DIMENSION B(78),C(12)
      C(1)=0.
      GO TO (1,2,3,4,5,6,7,8),I
    1 B(1)=DX
      RETURN
    2 C(2)=DX
      B(1)=DX
      B(2)=.5*DX
      B(3)=B(2)
      RETURN
    3 C(2)=.5*DX
      C(3)=DX
      B(1)=C(2)
      B(2)=-DX
      B(3)=DX+DX
      B(4)=.1666666667*DX
      B(5)=.6666666667*DX
      B(6)=B(4)
      RETURN
    4 C(2)=.01*DX
      C(3)=.6*DX
      C(4)=DX
      B(1)=C(2)
      B(2)=-17.46122448979*DX
      B(3)=18.06122448979*DX
      B(4)=59.69127516778*DX
      B(5)=-60.53065635308*DX
      B(6)=1.839381185303*DX
      B(7)=-2.5555555556*DX
      B(8)=2.853392683901*DX
      B(9)=.5767419962335*DX
      B(10)=.1254208754209*DX
      RETURN
    5 C(2)=.25*DX
      C(3)=C(2)
      C(4)=.5*DX
      C(5)=.75*DX
      C(6)=DX
      B(1)=C(2)
      B(2)=.125*DX
      B(3)=.125*DX
      B(4)=0.
      B(5)=-C(4)
      B(6)=C(6)
      B(7)=.1875*DX
      B(8)=0.
      B(9)=0.
      B(10)=.5625*DX
      B(11)=-.42857142857*DX
      B(12)=.285714285714*DX
      B(13)=1.7142857143*DX
      B(14)=-B(13)
      B(15)=1.14285714286*DX
      B(16)=.077777777778*DX
      B(17)=0.
      B(18)=.355555555556*DX
      B(19)=.13333333333*DX
      B(20)=B(18)
      B(21)=B(16)
      I=6
      RETURN
    6 C(2)=.1111111111*DX
      C(3)=.1666666667*DX
      C(4)=.3333333333*DX
      C(5)=.5*DX
      C(6)=.6666666667*DX
      C(7)=.8333333333*DX
      C(8)=DX
      B(1)=C(2)
      B(2)=.04166666667*DX
      B(3)=.125*DX
      B(4)=.1666666667*DX
      B(5)=-.5*DX
      B(6)=.6666666667*DX
      B(7)=-.625*DX
      B(8)=3.375*DX
      B(9)=-3.*DX
      B(10)=.75*DX
      B(11)=24.55555556*DX
      B(12)=-109.*DX
      B(13)=96.33333333*DX
      B(14)=-11.33333333*DX
      B(15)=.1111111111*DX
      B(16)=-3.8125*DX
      B(17)=14.125*DX
      B(18)=-9.833333333*DX
      B(19)=-1.375*DX
      B(20)=1.666666667*DX
      B(21)=.0625*DX
      B(22)=8.731707317*DX
      B(23)=-25.35365854*DX
      B(24)=12.21951219*DX
      B(25)=10.17073171*DX
      B(26)=-5.536585366*DX
      B(27)=-.1097560976*DX
      B(28)=.8780487805*DX
      B(29)=.04880952381*DX
      B(30)=0.
      B(31)=.2571428571*DX
      B(32)=.03214285714*DX
      B(33)=.3238095238*DX
      B(34)=B(32)
      B(35)=B(31)
      B(36)=B(29)
      I=8
      RETURN
    7 C(2)=.2222222222*DX
      C(3)=.3333333333*DX
      C(4)=.5*DX
      C(5)=.1666666667*DX
      C(6)=.8888888889*DX
      C(7)=.1111111111*DX
      C(8)=.8333333333*DX
      C(9)=DX
      B(1)=C(2)
      B(2)=.08333333333*DX
      B(3)=.25*DX
      B(4)=.125*DX
      B(5)=0.
      B(6)=.375*DX
      B(7)=.1064814814*DX
      B(8)=0.
      B(9)=.09722222222*DX
      B(10)=-.03703703704*DX
      B(11)=-5.673525377*DX
      B(12)=0.
      B(13)=-18.63374486*DX
      B(14)=7.22085048*DX
      B(15)=17.97530864*DX
      B(16)=.693329904*DX
      B(17)=0.
      B(18)=1.991769547*DX
      B(19)=-.7105624143*DX
      B(20)=-1.874643875*DX
      B(21)=.01121794872*DX
      B(22)=-.5634259259*DX
      B(23)=0.
      B(24)=-2.013888889*DX
      B(25)=1.261073318*DX
      B(26)=1.851282051*DX
      B(27)=.05951726845*DX
      B(28)=.2387755102*DX
      B(29)=.09356936416*DX
      B(30)=0.
      B(31)=-.4855491329*DX
      B(32)=-.08092485549*DX
      B(33)=2.761227212*DX
      B(34)=-.3964976497*DX
      B(35)=-1.852251794*DX
      B(36)=.9604268564*DX
      B(37)=.05148809524*DX
      B(38)=0.
      B(39)=0.
      B(40)=.3587949466*DX
      B(41)=.2967032967*DX
      B(42)=-.02758886522*DX
      B(43)=B(42)
      B(44)=B(41)
      B(45)=B(37)
      I=9
      RETURN
    8 C( 2)= .11111111111*DX
      C( 3)= .16666666667*DX
      C( 4)= .25*DX
      C( 5)= .1*DX
      C( 6)= C(3)
      C( 7)= .5*DX
      C( 8)= .666666666667*DX
      C( 9)= .33333333333*DX
      C(10)= .83333333333*DX
      C(11)= C(10)
      C(12)= DX
      B( 1)= C(2)
      B( 2)= .041666666667*DX
      B( 3)= .125*DX
      B( 4)= .0625*DX
      B( 5)= 0.
      B( 6)= .1875*DX
      B( 7)= .058*DX
      B( 8)= 0.
      B( 9)= .066*DX
      B(10)= -.024*DX
      B(11)= .033950617284*DX
      B(12)= 0.
      B(13)= 0.
      B(14)= .0041152263374*DX
      B(15)= .12860082305*DX
      B(16)= -.58333333333*DX
      B(17)= 0.
      B(18)= 0.
      B(19)= 2.1111111111*DX
      B(20)= 3.4722222222*DX
      B(21)= -4.5*DX
      B(22)= -.12345678901*DX
      B(23)= 0.
      B(24)= 0.
      B(25)= -.1316872428*DX
      B(26)= .51440329218*DX
      B(27)= 0.
      B(28)= .40740740741*DX
      B(29)= 3.6265432099*DX
      B(30)= 0.
      B(31)= 0.
      B(32)= -10.666666667*DX
      B(33)= -19.290123457*DX
      B(34)= 26.*DX
      B(35)= .74691358025*DX
      B(36)= -.083333333333*DX
      B(37)= .90432098765*DX
      B(38)= 0.
      B(39)= 0.
      B(40)= -2.6296296296*DX
      B(41)= -4.2438271605*DX
      B(42)= 5.6666666667*DX
      B(43)= -.36419753086*DX
      B(44)= .5*DX
      B(45)= DX
      B(46)= .80432098765*DX
      B(47)= 0.
      B(48)= 0.
      B(49)= -2.6296296296*DX
      B(50)= -4.2438271605*DX
      B(51)= 6.1666666667*DX
      B(52)= .63580246914*DX
      B(53)= 0.
      B(54)= 0.
      B(55)= .1*DX
      B(56)= -1.9410569106*DX
      B(57)= 0.
      B(58)= 0.
      B(59)= 6.9376693767*DX
      B(60)= 11.009485095*DX
      B(61)= -14.926829268*DX
      B(62)= .085365853659*DX
      B(63)= -.16463414634*DX
      B(64)= -.43902439024*DX
      B(65)= -.29268292683*DX
      B(66)= .73170731707*DX
      B(67)= .04880952381*DX
      B(68)= 0.
      B(69)= 0.
      B(70)= 0.
      B(71)= 0.
      B(72)= .25714285714*DX
      B(73)= .32380952381*DX
      B(74)= .032142857143*DX
      B(75)= B(74)
      B(76)= .042857142857*DX
      B(77)= .21428571429*DX
      B(78)=B(67)
      I=12
      RETURN
      END
      SUBROUTINE STEPS(EPS,VERTNO)
C     UNALTERED DPFG VERSION
      implicit complex(a-h,o-z)
      COMMON/C/STEP
      DIMENSION STEP(8)
      PS=ALOG(cabs(EPS))
      VERTNO=-PS
      FAC=1.
      DO 2 N=1,8
      FN=N+1
      FAC=FAC*FN
      X=(ALOG(cabs(FAC))+PS)/FN
      X=EXP(X)
      S=X
      DO 1 I=1,N
    1 S=X*EXP(-S/FN)
    2 STEP(N)=S
      RETURN
      END
      SUBROUTINE MATADJ(C,CA)
      implicit complex(a-h,o-z)
      DIMENSION C(3,3),CA(3,3)
      CA(1,1)=C(2,2)*C(3,3)-C(3,2)*C(2,3)
      CA(1,2)=C(3,2)*C(1,3)-C(1,2)*C(3,3)
      CA(1,3)=C(1,2)*C(2,3)-C(2,2)*C(1,3)
      CA(2,1)=C(3,1)*C(2,3)-C(2,1)*C(3,3)
      CA(2,2)=C(1,1)*C(3,3)-C(3,1)*C(1,3)
      CA(2,3)=C(2,1)*C(1,3)-C(1,1)*C(2,3)
      CA(3,1)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
      CA(3,2)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
      CA(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
      RETURN
      END
      SUBROUTINE SPSJFG(I,A)
C     UNALTERED DPFG VERSION.
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/J/KG
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),A(6,3)
      X=R(I)
      RO=RHO(I)
      GR=G(I)
      P=VP(I)
      VPSQ=P*P
      P=VS(I)
      VSSQ=P*P
      ZETA=4.*CON*RO
      XI=GR/X
      XSQ=X*X
      ALFSQ=(WSQ+ZETA+XI)/VPSQ
      BETASQ=WSQ/VSSQ
      GAMSQ=4.*FL3*XI*XI/(VPSQ*VSSQ)
      DELSQ=SQRT((BETASQ-ALFSQ)*(BETASQ-ALFSQ)+GAMSQ)
      FKSQ=.5*(ALFSQ+BETASQ+DELSQ)
      QSQ=FKSQ-DELSQ
      XL=X**L
      XLP1=XL*X
      XLP2=XLP1*X
      XLM1=XL/X
      FKXSQ=FKSQ*XSQ
      QXSQ=QSQ*XSQ
      QK=QXSQ*FKXSQ
      QPK=QXSQ+FKXSQ
      FLU=RO*VPSQ
      FU=RO*VSSQ
      D0=1.
      H0=-VPSQ/(FL1*VSSQ)
      K=1
    5 C=2.
      B=FL2+2.
      F=2.
      C2=1./(C*B)
      D1=C2*(FL3*XI*H0/VPSQ-ALFSQ*D0)*XSQ
      H1=C2*(XI*D0/VSSQ-BETASQ*H0)*XSQ
      U=C2*(FL3*H0+(FL+F)*D0)
      V=C2*(D0+(FL1+F)*H0)
      P=C2*D0
      S=(FL2+F)*P-U
      H=H0+H1
      D=D0+D1
    6 C1=C2
      C=C+2.
      B=B+2.
      F=F+2.
      C2=1./(C*B)
      UN=C2*(FL3*H1+(FL+F)*D1)
      VN=C2*(D1+(FL1+F)*H1)
      PN=C2*D1
      SN=(FL2+F)*PN-UN
      D2=-C2*(D1*ALFSQ-H1*FL3*XI/VPSQ)*XSQ
      H2=-C2*(H1*BETASQ-D1*XI/VSSQ)*XSQ
      D=D+D2
      H=H+H2
      D0=D1
      D1=D2
      H0=H1
      H1=H2
      U=U+UN
      V=V+VN
      P=P+PN
      S=S+SN
      TE=ABS (D2/D)
      IF(real(TE)-real(EPS)) 7,6,6
    7 TE=ABS(H2/H)
      IF(real(TE)-real(EPS)) 8,6,6
    8 C=C+2.
      B=B+2.
      F=F+2.
      C2=1./(C*B)
      UN=C2*(FL3*H1+(FL+F)*D1)
      VN=C2*(D1+(FL1+F)*H1)
      PN=C2*D1
      SN=(FL2+F)*PN-UN
      A(1,K)=(U+UN)*XLP1
      A(2,K)=(V+VN)*XLP1
      A(3,K)=FLU*D*XL+2.*FU*(FL3*A(2,K)-2.*A(1,K))/X
      A(4,K)=FU*(H*XL+2.*(A(1,K)-A(2,K))/X)
      A(5,K)=ZETA*(P+PN)*XLP2
      A(6,K)=ZETA*(S+SN)*XLP1
      GO TO (9,10),K
    9 K=2
      D0=0.
      H0=-1.
      GO TO 5
   10 A(1,3)=XLM1*FL
      A(2,3)=XLM1
      A(3,3)=2.*FU*(FL3*A(2,3)-2.*A(1,3))/X
      A(4,3)=2.*FU*(A(1,3)-A(2,3))/X
      B=XI*FL-WSQ
      A(5,3)=B*XL
      A(6,3)=(FL2*B-ZETA*FL)*XLM1
      A(5,2)=A(5,2)+FL1*VSSQ*XL
      A(6,2)=A(6,2)+FL2*FL1*VSSQ*XLM1
      IS=I
      DO 20 I=1,3,2
      A(2,I)=SFL3*A(2,I)
   20 A(4,I)=SFL3*A(4,I)
      DO 21 I=1,5,2
   21 A(I,2)=A(I,2)/SFL3
      A(6,2)=A(6,2)/SFL3
      I=IS
      IF(KG.EQ.1) RETURN
      DO 30 J=1,6
   30 A(J,3)=0.
      A(5,1)=0.
      A(5,2)=0.
      A(6,1)=0.
      A(6,2)=0.
      RETURN
      END

      SUBROUTINE CCOEF(X,I,A)
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),A(2,2)
      J=I
    1 J=J-1
      IF(real(R(J))-real(R(I))) 2,1,1
    2 DEL=(X-R(J))/(R(I)-R(J))
      RO=RHO(J)+DEL*(RHO(I)-RHO(J))
      GR=G(J)+DEL*(G(I)-G(J))
      ROPR=(RHO(I)-RHO(J))/(R(I)-R(J))
      Z=1./X
      ZSQ=Z*Z
      A(1,1)=0.
      A(1,2)=1.
      A(2,1)=FL3*ZSQ+4.*CON*ROPR/GR
      A(2,2)=-2.*Z
      RETURN
      END

      SUBROUTINE CMTR(X,Y,U,V)
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/C/STEP
      COMMON/I/B,C,DX,IN
      COMMON/E/NDIS,KDIS
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),F(2),A(2,2),S(2),C(12),B(78),H(12,2),STEP(8),KDIS(28)
      IF(real(Y)-real(X)) 1,1,2
    1 RETURN
    2 YS=Y
      XS=X
      NOC=KDIS(2)
      F(1)=U
      F(2)=V
      I=2
    3 IF(real(Y)-real(R(I))) 5,5,4
    4 I=I+1
      IF(I-NOC) 3,5,5
    5 RO=RHO(I)
      FLAMB=FLAM(I)
    6 QSQ=ABS(4.*CON*RO*RO/FLAMB-FL3/(X*X))
      Q=SQRT(QSQ)+1./X
      DX=STEP(8)/Q
      Y=X+DX
      IF(real(Y)-real(YS)) 11,11,10
   10 Y=YS
   11 DX=Y-X
      DS=Q*DX
      DO 13 J=1,7
      IF(real(DS).LE.real(STEP(J))) GO TO 12
      GO TO 13
   12 IN=J
      GO TO 14
   13 CONTINUE
      IN=8
   14 CALL SHANKS
      S(1)=F(1)
      S(2)=F(2)
      DO 25 NI=1,IN
      Z=X+C(NI)
      CALL CCOEF(Z,I,A)
      H(NI,1)=A(1,1)*F(1)+A(1,2)*F(2)
      H(NI,2)=A(2,1)*F(1)+A(2,2)*F(2)
      F(1)=S(1)
      F(2)=S(2)
      DO 25 M=1,NI
      K1=M+NI*(NI-1)/2
      F(1)=F(1)+B(K1)*H(M,1)
   25 F(2)=F(2)+B(K1)*H(M,2)
      X=Y
      IF(real(Y)-real(YS)) 6,26,26
   26 X=XS
      Y=YS
      U=F(1)
      V=F(2)
      RETURN
      END
      SUBROUTINE SCOEF(X,I,C)
      parameter(numb=1000)
      implicit complex(a-h,o-z)
C     CONVENTION USED IS DELSQ PHI = 4 PI G RHO.
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),A(6,6),B(6,6),C(6,6)
      J=I
    1 J=J-1
      IF(real(R(J))-real(R(I))) 2,1,1
    2 DEL=(X-R(J))/(R(I)-R(J))
      RO=RHO(J)+DEL*(RHO(I)-RHO(J))
      FU=FMU(J)+DEL*(FMU(I)-FMU(J))
      FLU=FLAM(J)+DEL*(FLAM(I)-FLAM(J))
      GR=G(J)+DEL*(G(I)-G(J))
      DO 3 K=1,6
      DO 3 J=1,6
      A(J,K)=0.
      B(J,K)=0.
    3 C(J,K)=0.
      D=1./(FLU+FU+FU)
      E=(3.*FLU+FU+FU)*D*FU
      ZETA=4.*CON*RO
      A(3,1)=4.*E
      A(4,1)=-2.*E
      A(3,2)=A(4,1)*FL3
      A(4,2)=2.*FU*(2.*(FLU+FU)*FL3*D-1.)
      B(1,1)=-2.*FLU*D
      B(2,1)=-1.
      B(3,1)=-4.*RO*GR
      B(4,1)=RO*GR
      B(6,1)=ZETA*FL1
      B(1,2)=FLU*D*FL3
      B(2,2)=1.
      B(3,2)=B(4,1)*FL3
      B(6,2)=-ZETA*FL3
      B(3,3)=-4.*FU*D
      B(4,3)=-FLU*D
      B(3,4)=FL3
      B(4,4)=-3.
      B(3,5)=RO*FL1
      B(4,5)=-RO
      B(5,5)=-FL1
      B(6,6)=FL-1.
      C(3,1)=-RO*WSQ
      C(5,1)=ZETA
      C(4,2)=C(3,1)
      C(1,3)=D
      C(2,4)=1./FU
      C(3,6)=-RO
      C(5,6)=1.
      Z=1./X
      ZSQ=Z*Z
      C(1,1)=Z*B(1,1)
      C(2,1)=Z*B(2,1)
      C(3,1)=ZSQ*A(3,1)+Z*B(3,1)+C(3,1)
      C(4,1)=ZSQ*A(4,1)+Z*B(4,1)
      C(6,1)=Z*B(6,1)
      C(1,2)=Z*B(1,2)
      C(2,2)=Z*B(2,2)
      C(3,2)=ZSQ*A(3,2)+Z*B(3,2)
      C(4,2)=ZSQ*A(4,2)+C(4,2)
      C(6,2)=Z*B(6,2)
      C(3,3)=Z*B(3,3)
      C(4,3)=Z*B(4,3)
      C(3,4)=Z*B(3,4)
      C(4,4)=Z*B(4,4)
      C(3,5)=Z*B(3,5)
      C(4,5)=Z*B(4,5)
      C(5,5)=Z*B(5,5)
      C(6,6)=Z*B(6,6)
      C(3,5)=-C(3,5)
      C(3,6)=-C(3,6)
      C(4,5)=-C(4,5)
      C(5,1)=-C(5,1)
      C(6,1)=-C(6,1)
      C(6,2)=-C(6,2)
      RETURN
      END

      SUBROUTINE SMTR(X,Y,F)
      parameter(numb=1000)
      implicit complex(a-h,o-z)
      COMMON CON,R,RHO,VP,VS,VFSQ,FLAM,FMU,G,RHOBAR,W,WSQ,
     1DEL1,DEL2,QS,QF,Q,EPS,SFL3,FORD,VERTNO,PI,FL,
     2FL1,FL2,FL3,N1,N2,N,JCOM,L
      COMMON/C/STEP
      COMMON/I/B,C,DX,IN
      COMMON/J/KG
      DIMENSION R(numb),RHO(numb),VP(numb),VS(numb),VFSQ(numb),
     2  FLAM(numb),FMU(numb),
     1 G(numb),F(6,3),A(6,6),S(6,3),C(12),B(78),H(12,6,3),STEP(8)
      IF(real(Y)-real(X)) 1,1,2
    1 RETURN
    2 YS=Y
      XS=X
      KK=4-KG
      JJ=12+KG*(2*KG-8)
      I=2
    3 IF(real(Y)-real(R(I))) 5,5,4
    4 I=I+1
      IF(I-N) 3,5,5
    5 V=VP(I)
      VPSQ=V*V
      VSSQ=VS(I)*VS(I)
      ZETA=4.*CON*RHO(I)
      RI=R(I)
      XI=G(I)/RI
      ALFSQ=(WSQ+ZETA+XI)/VPSQ
      BETASQ=WSQ/VSSQ
      GAMSQ=4.*FL3*XI*XI/(VSSQ*VPSQ)
      DELSQ=SQRT((BETASQ-ALFSQ)*(BETASQ-ALFSQ)+GAMSQ)
      FKSQ=.5*(ALFSQ+BETASQ+DELSQ)
      QSQ=FKSQ-DELSQ
      SFL3=SQRT(FL3)
    6 Q=SFL3/X
      QS=SQRT(ABS(FKSQ-FL3/(X*X)))+1./X
      QF=SQRT(ABS(QSQ-FL3/(X*X)))+1./X
      GO TO (61,7,60),KG
   60 Q=QS+QF
      GO TO 9
   61 IF(real(Q)-real(QF)) 7,8,8
    7 Q=QF
    8 IF(real(Q)-real(QS)) 81,9,9
   81 Q=QS
    9 DX=STEP(8)/Q
      Y=X+DX
      IF(real(Y)-real(YS)) 11,11,10
   10 Y=YS
   11 DX=Y-X
      DS=Q*DX
      DO 13 J=1,7
      IF(real(DS).LE.real(STEP(J))) GO TO 12
      GO TO 13
   12 IN=J
      GO TO 14
   13 CONTINUE
      IN=8
   14 CALL SHANKS
      DO 25 J=1,JJ
      DO 25 K=1,KK
   25 S(J,K)=F(J,K)
      DO 27 NI=1,IN
      Z=X+C(NI)
      CALL SCOEF(Z,I,A)
      DO 26 J=1,JJ
      DO 26 K=1,KK
      H(NI,J,K)=0.
      DO 26 M=1,JJ
   26 H(NI,J,K)=H(NI,J,K)+A(J,M)*F(M,K)
      DO 27 J=1,JJ
      DO 27 K=1,KK
      F(J,K)=S(J,K)
      DO 27 M=1,NI
      K1=M+NI*(NI-1)/2
   27 F(J,K)=F(J,K)+B(K1)*H(M,J,K)
      X=Y
      IF(real(Y)-real(YS)) 6,29,29
   29 X=XS
      Y=YS
      RETURN
      END

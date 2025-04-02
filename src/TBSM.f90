!----------------------------------------------------------------
! CTEMsoft                                                      !
! Copyright (c) 2001, Marc De Graef/Carnegie Mellon University  !
! All rights reserved.                                          !
!                                                               !
! The source code in this file is protected by a BSD license;   !
! for the complete text of the license see the file BSD.license.!
!                                                               !
!----------------------------------------------------------------
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "TBSM.f90"
!                                    created: 4/28/01  {9:29:46 AM} 
!                                last update: 4/28/01 {10:37:32 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: Two beam scattering matrix program 
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  4/28/01  MDG 1.0 original
!  5/27/01  MDG 2.0 f90
! ###################################################################
program TBSM

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

integer     :: n(1),g(3),ns,nt
real        :: x(1),wmax

 progname = 'TBSM.f90'
 progdesc = 'Two-beam computations: comparing scattering matrix and analytical'
 call CTEMsoft
                                                                                                  

! first get the crystal data and microscope voltage
 SG % SYM_reduce=.TRUE.
 call CrystalData
 call GetVoltage
! generate all atom positions
 call CalcPositions('v')
! get the reciprocal lattice vector
 mess = 'Diffraction vector :'; call Message("(A)")
 call GetIndex(g,'r')
! some more parameters
 mess = 'Enter maximum value of w: '; call GetReal(1)
 wmax = io_real(1)
 mess = 'Maximum foil thickness : '; call GetReal(1)
! do the computation
 call CalcTBSM(g,wmax,io_real(1))
end program



subroutine CalcTBSM(g,wmax,tmax)

use local
use io
use crystal
use crystalvars
use diffraction
use constants
use postscript
use dynamical

integer,parameter      :: ns=512, nt=512
integer(kind=irg)      :: g(3),ind(3)
real(kind=sgl)         :: qr,qi,bg,Vmod,Vphase,Vpmod,Vpphase,wmax,xg,xgp,xgpz, &
                          tmax, Ar(2,2), Ai(2,2), sg, dz, dsg, p(2),q(2), &
                          BF(512,512,2),DF(512,512,2),It,Is

! normal aborption factor
 ind = [0,0,0]
 call CalcUcg(ind)
 xgpz= rlp%xgp
! extinction distance and anomalous absorption length
 call CalcUcg(g)
 xgp = rlp%xgp
 xg  = rlp%xg
 imanum=0
 dz = tmax/float(nt-1)
! compute the scattering matrix SM for each orientation
 dsg = 2.0*wmax/float(ns-1)
 mess = 'Starting scattering matrix computation'; call Message("(A)")
 mess = 'Number of computations per * = '; oi_int(1)=8*ns; call WriteInt(1,"(I4)")
! allocate arrays
 allocate(SMr(512,2,2))
 allocate(SMi(512,2,2))
 allocate(phir(512,2))
 allocate(phii(512,2))
 T0 = SECNDS(0.0)
 do i=1,ns
  sg = (-wmax+dsg*float(i))/xg
  call TBCalcSM(Ar,Ai,sg,dz,xg,xgp,xgpz,bg)
  SMr(i,1:2,1:2)=Ar(1:2,1:2)
  SMi(i,1:2,1:2)=Ai(1:2,1:2)
 end do
! compute intensities for the first row
  do i=1,ns
   BF(i,1,1)=SMr(i,1,1)**2+SMi(i,1,1)**2
   DF(i,1,1)=SMr(i,2,1)**2+SMi(i,2,1)**2
! the first column of SM is the actual initial wavefunction at z=dz
   phir(i,1)=SMr(i,1,1)
   phir(i,2)=SMr(i,2,1)
   phii(i,1)=SMi(i,1,1)
   phii(i,2)=SMi(i,2,1)
  end do
! loop to maximum thickness
  do l=2,nt
   do k=1,ns
    do i=1,2
     p(i)=0.D0
     q(i)=0.D0
     do j=1,2
      p(i)=p(i)+SMr(k,i,j)*phir(k,j)-SMi(k,i,j)*phii(k,j)
      q(i)=q(i)+SMi(k,i,j)*phir(k,j)+SMr(k,i,j)*phii(k,j)
     end do
    end do
    phir(k,1:2)=p(1:2)
    phii(k,1:2)=q(1:2)
   end do
   if (mod(l,8).eq.0) then
     mess = '*'
     call Message("(1A,$)")
   end if
   do j=1,ns
    BF(j,l,1) = phir(j,1)**2+phii(j,1)**2
    DF(j,l,1) = phir(j,2)**2+phii(j,2)**2
   end do
  end do
  T1 = SECNDS(T0)
  mess = 'done'; call Message("(A)")
  mess = 'Total computation time [s] '; oi_real(1)=T1; call WriteReal(1,"(F)")
  mess = 'Starting direct analytical solution'; call Message("(/A)")
! next, redo the computation, but this time use TBCalcInten
  T0 = SECNDS(0.0)
  do i=1,ns
   sg = (-wmax+dsg*float(i))/xg
   do j=1,nt
    t = float(j)*dz
    call TBCalcInten(It,Is,sg,t,xg,xgp,xgpz,bg)
    BF(i,j,2) = It
    DF(i,j,2) = Is
   end do
   if (mod(i,8).eq.0) then
     mess = '*'
     call Message("(1A,$)")
   end if
  end do
  T2 = SECNDS(T0)
  mess = 'done'; call Message("(A)")
  mess = 'Total computation time [s] '; oi_real(1)=T2; call WriteReal(1,"(F)")
! compare the two computations
  bfdiff = 0.0
  dfdiff = 0.0
  do i=1,ns
   do j=1,nt
    bfdiff = bfdiff + (BF(i,j,1)-BF(i,j,2))**2
    dfdiff = dfdiff + (DF(i,j,1)-DF(i,j,2))**2
   end do
  end do
  bfdiff = bfdiff/float(ns)/float(nt)
  dfdiff = dfdiff/float(ns)/float(nt)
  mess = 'Average difference in BF '; oi_real(1)=bfdiff; call WriteReal(1,"(F)")
  mess = 'Average difference in DF '; oi_real(1)=dfdiff; call WriteReal(1,"(F)")
  mess = 'Images computed, preparing for PS output'; call Message("(A)")
! open PostScript file
  call PS_openfile
  call PS_newpage(.TRUE.,'Two Beam Rocking Curves')
  call PS_setfont(PSfonts(2),0.16)
  call PS_text(2.5,8.6,'Scattering Matrix Solution')
  call PS_text(2.5,2.1,'Direct Analytical Solution')
  call PS_setfont(PSfonts(2),0.10)
  call PS_text(0.5,1.8,'Input file')
  call PS_text(2.6,1.8,cell % fname)
  call PS_text(0.5,1.6,'Active reflection')
  call PS_textint(2.5,1.6,' ',g(1))
  call PS_textint(2.6,1.6,' ',g(2))
  call PS_textint(2.7,1.6,' ',g(3))
  call PS_text(0.5,1.4,'Extinction distance  [nm]')
  call PS_textvar(2.5,1.4,' ',xg)
  call PS_text(0.5,1.2,'Anomalous absorption length [nm]')
  call PS_textvar(2.5,1.2,' ',xgp)
  call PS_text(0.5,1.0,'Normal absorption length [nm]')
  call PS_textvar(2.5,1.0,' ',xgpz)
  call PS_text(0.5,0.8,'Accelerating Voltage [V]')
  call PS_textvar(2.5,0.8,' ',sngl(mAccvol))
  call PS_text(0.5,0.6,'Maximum w')
  call PS_textvar(2.5,0.6,' ',wmax)
  call PS_text(0.5,0.4,'Maximum thickness [nm]')
  call PS_textvar(2.5,0.4,' ',tmax)
! Bright Field image
  allocate(imaint(512,512))
  imaint(1:ns,1:nt)=int(255.0*BF(1:ns,1:nt,1))
  x0=0.2
  y0=5.5
  npx=ns
  npy=nt
  scl=3.0
  call PS_DumpImage(x0,y0,npx,npy,scl)
! Dark Field image
  imaint(1:ns,1:nt)=int(255.0*DF(1:ns,1:nt,1))
  x0=3.4
  y0=5.5
  npx=ns
  npy=nt
  scl=3.0
  call PS_DumpImage(x0,y0,npx,npy,scl)
! Bright Field image
  imaint(1:ns,1:nt)=int(255.0*BF(1:ns,1:nt,2))
  x0=0.2
  y0=2.3
  npx=ns
  npy=nt
  scl=3.0
  call PS_DumpImage(x0,y0,npx,npy,scl)
! Dark Field image
  imaint(1:ns,1:nt)=int(255.0*DF(1:ns,1:nt,2))
  x0=3.4
  y0=2.3
  npx=ns
  npy=nt
  scl=3.0
  call PS_DumpImage(x0,y0,npx,npy,scl)
! close Postscript file
  call PS_closefile
end subroutine
       


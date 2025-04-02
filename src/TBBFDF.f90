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
!  FILE: "TBBFDF.f90"
!                                    created: 12/11/98 {9:29:46 AM} 
!                                last update: 4/18/01 {10:37:32 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: Two beam bright field dark field image pairs; the
!               file TBBFDF.routines contains a number of options
!               which can be pasted into this program.  See instructions
!               in BFDF.routines.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  12/11/98 MDG 1.0 original
!   5/27/01 MDG 2.0 f90
! ###################################################################
program tbBFDF

use local
use io
use constants
use files
use crystal
use symmetry
use diffraction
use postscript
use dynamical

real,parameter      :: thr = 1.0E-6
real,allocatable    :: BF(:,:), DF(:,:), exe(:,:), z(:,:)
real                :: dz, xirat ,zmin, zmax
integer             :: jdim,zero(3),ind(3),seed
real                :: qr,qi,bg,r(200),p(200),xig,xigp,betag,xizero

 progname = 'TBBFDF.f90'
 progdesc = 'Two-beam bright field-dark field images, using BFDF.routines'
 call CTEMsoft
                                                                               
! compute the relevant parameters for a given crystal
 imanum = 0
 call CrystalData
 call GetVoltage
 call CalcPositions('v')
! initialize all arrays
 mess = 'Enter dimension of image array (<= 512) '; call GetInt(1)
 jdim = min(io_int(1),512)
 allocate(BF(jdim,jdim))
 allocate(DF(jdim,jdim))
 allocate(exe(jdim,jdim))
 allocate(z(jdim,jdim))
! the following array was declared in module postscript.f90
 allocate(imaint(jdim,jdim))
! extinction distance and absorption length for given plane
 mess = 'Indices for reciprocal lattice vector g'; call Message("(A)")
 call GetIndex(ind,'r')
 call CalcUcg(ind)
 xigp = rlp%xgp
 xig  = rlp%xg
 betag = rlp%Vpphase-rlp%Vphase
 mess = 'Extinction distance [nm] = '; oi_real(1) = xig; call WriteReal(1,"(f10.4)")
 mess = 'Absorption length   [nm] = '; oi_real(1) = xigp; call WriteReal(1,"(f10.4)")
! normal absorption length
 zero(1)=0
 zero(2)=0
 zero(3)=0
 call CalcUcg(zero)
 xizero = rlp%xgp
 mess = 'Normal Absorption   [nm] = '; oi_real = sngl(xizero); call WriteReal(1,"(f10.4)")
! the following lines are taken from the BFDF.routines file
! cut from ! REPLACE START until the line   ! REPLACE STOP
! 
! REPLACE START
 T0 = SECNDS(0.0)
 mess = 'Average thickness, minimum thickness [nm]  = '; call Message("(A)")
 mess = '(minimum thickness may be <0, to create holes)'; call GetReal(2)
 zav = io_real(1)
 zdev = io_real(2)
 mess = 'How many random waves ? (<100) = '; call GetInt(1)
 nw = oi_int(1)
 nw = min(nw,100)
 ff = cPi/float(jdim)
! superimpose a bunch of waves with random periodicities and phases
 T1 = SECNDS(T0)
 seed = int(1000.0*T1)
 do k=1,2*nw
  r(k) = ran(seed)*10.0
  p(k) = ran(seed)*3.0
 end do
 do i=1,jdim
  fi = float(i)*ff
  do j=1,jdim
   fj=float(j)*ff
   z(i,j) = 0.0 
   do k=1,nw
    z(i,j) = z(i,j) + cos(r(2*k)*fi+p(2*k))*cos(r(2*k-1)*(fi+fj)+p(2*k-1)) + &
                      sin(r(2*k-1)*fj+p(2*k))*sin(r(2*k)*(fi-fj)+p(2*k-1))
   end do
  end do
 end do
 zmin = minval(z)
 zmax = maxval(z)
 ztot = sum(z)/float(jdim)**2
 oi_real(1) = zmin
 oi_real(2) = zmax
 mess = 'min max after generation : '; call WriteReal(2,"(1x,2(f8.3,2x))")
! create the excitation error array 
 q = (zav-zdev)/(ztot-zmin)
 do i=1,jdim
  ss = (-0.5+1.0*float(i-1)/float(jdim))/xig
  do j=1,jdim
   exe(i,j) = ss
   z(i,j) = zdev+q*(z(i,j)-zmin)
   if (z(i,j).lt.0.0) z(i,j)=0.0
  end do
 end do
 zmi = minval(z)
 zma = maxval(z)
 ztot = sum(z)/float(jdim)**2
 oi_real(1) = zmi
 oi_real(2) = zma
 mess = 'min max after scaling    : '; call WriteReal(2,"(1x,2(f8.3,2x))")
 mess = 'Actual average thickness = '; oi_real(1)=ztot; call WriteReal(1,"(f10.4)")
 zmin=zmi 
 zmax=zma 
! REPLACE STOP
!
! compute intensities
 T0= SECNDS(0.0)
 mess = 'computation start '; call Message("(A)")
 iflag=0
 do i=1,jdim
  do j=1,jdim
   call TBCalcInten(BF(i,j),DF(i,j),exe(i,j),z(i,j),xig,xigp,xizero,betag)
   if (iflag.eq.0) then 
    if ((BF(i,j).gt.1.D0).or.(DF(i,j).gt.1.D0)) iflag=1
   endif
  end do
 end do
 T1= SECNDS(T0)
 mess = 'computation end   -> '; oi_real(1)=T1; call WriteReal(1,"(F,' seconds')")
 tps = T1/float(jdim)**2
 mess = 'time per pixel    -> '; oi_real(1)=tps; call WriteReal(1,"(F,' seconds')")
 if (iflag.eq.1) then 
  mess = 'WARNING: one or more intensities larger than 1 !'; call Message("(A)")
  mess = 'Normal absorption length is too large'; CALL MESSAGE("(A)")
 endif
! open PostScript file
 call PS_openfile
 pspage = 0
 call PS_newpage(.TRUE.,'Two Beam Bright Field - Dark Field')
 call PS_setfont(PSfonts(2),0.10)
 call PS_text(0.5,1.8,'Input file')
 call PS_text(2.6,1.8,PS % psname)
 call PS_text(0.5,1.6,'Active reflection')
 call PS_textint(2.5,1.6,' ',ind(1))
 call PS_textint(2.6,1.6,' ',ind(2))
 call PS_textint(2.7,1.6,' ',ind(3))
 call PS_text(0.5,1.4,'Extinction distance  [nm]')
 call PS_textvar(2.5,1.4,' ',xig)
 call PS_text(0.5,1.2,'Anomalous absorption length [nm]')
 call PS_textvar(2.5,1.2,' ',xigp)
 call PS_text(0.5,1.0,'Normal absorption length [nm]')
 call PS_textvar(2.5,1.0,' ',xizero)
 call PS_text(0.5,0.8,'Accelerating Voltage [V]')
 call PS_textvar(2.5,0.8,' ',sngl(mAccvol))
 call PS_text(0.5,0.6,'Minimum thickness [nm]')
 call PS_textvar(2.5,0.6,' ',zmin)
 call PS_text(0.5,0.4,'Maximum thickness [nm]')
 call PS_textvar(2.5,0.4,' ',zmax)
! Bright Field image
 do i=1,jdim
  do j=1,jdim
   imaint(i,j)=int(255.0*BF(i,j))
  end do
 end do 
 x0=0.2
 y0=4.0
 npx=jdim
 npy=jdim
 scl=3.0
 call PS_DumpImage(x0,y0,npx,npy,scl)
! Dark Field image
 do i=1,jdim
  do j=1,jdim
   imaint(i,j)=int(255.0*DF(i,j))
  end do
 end do
 x0=3.4
 y0=4.0
 npx=jdim
 npy=jdim
 scl=3.0
 call PS_DumpImage(x0,y0,npx,npy,scl)
! close Postscript file
 call PS_closefile
end program

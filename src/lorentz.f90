!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:lorentz.f90                                                          !
! Copyright (c) 2001, 2002  Marc De Graef/Carnegie Mellon University (CMU)      !
!                                                                               !
!     This program is free software; you can redistribute it and/or modify      !
!     it under the terms of the GNU General Public License as published by      !
!     the Free Software Foundation; either version 2 of the License, or         !
!     (at your option) any later version.                                       !
!                                                                               !
!     This program is distributed in the hope that it will be useful,           !
!     but WITHOUT ANY WARRANTY; without even the implied warranty of            !
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
!     GNU General Public License for more details.                              !
!                                                                               !
!     You should have received a copy of the GNU General Public License         !
!     along with this program; if not, write to the Free Software               !
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA !
!                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "lorentz.f90"
!                                    created:  4/16/97 
! 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description:  This module contains routines used for Lorentz microscopy
!                image simulations.
!                The phasemap routine uses the equations derived by M. Mansuripur
!                J. Appl. Phys., vol. 69, pp. 2455-2464 (1991)
!                to compute the Aharonov-Bohm phase shift caused
!                by a 2-D magnetisation distribution. 
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   4/16/97 MDG 1.0 original (Fortran-77)
!   9/21/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################

module lorentz_mod

use local

IMPLICIT NONE

complex,allocatable      :: LCTF(:,:)


contains

subroutine PhaseMap(jdim,B0t,beam,fname,gname)

use local
use constants
use files
use io
use TIFF_global
use TIFF_f90

IMPLICIT NONE

! parameters for fftw package
integer(kind=irg),parameter         ::  FFTW_FORWARD=-1,FFTW_BACKWARD=1,  &
                                        FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1, &
                                        FFTW_ESTIMATE=0,FFTW_MEASURE=1, &
                                        FFTW_OUT_OF_PLACE=0,FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16, &
                                        FFTW_THREADSAFE=128

character(20)             :: fname,gname
integer(kind=irg)         :: i,j,jdim,jdim2,jdim4,nn(2),ix,iy,iloc,izero,plan
logical                   :: igp
real(kind=sgl)            :: fjdim,prefac,d1r,d1i,d2r,d2i,d3r,d3i,d4r,d4i,beam(3),b,sx,sy,s,sigx,sigy,gp, &
                             psig,pre,B0t,mi,ma
real(kind=sgl),allocatable:: bx(:),by(:),bz(:),mx(:)

INTENT(IN)                :: jdim,B0t
INTENT(INOUT)             :: beam,fname
          
  fjdim=1.0/float(jdim)
  jdim2 = 2*jdim
  jdim4 = 2*jdim**2
  prefac= 2.0*cCharge*1.0E-18/cPlanck*fjdim**2
  allocate(bx(jdim4),by(jdim4),bz(jdim4),mx(jdim**2))

! normalize the incident beam direction
  b = sqrt(sum(beam**2))
  beam=beam/b

write (*,*) 'Normalized beam direction ',beam

! if the beam is along the z-axis, we do not need to compute the
! function G_p(ts) (controlled by variable igp)
  igp=.TRUE.
  if ((beam(1).eq.0.0).AND.(beam(2).eq.0.0)) then
   igp=.FALSE.
  endif

! read the magnetization pattern from file
  mess = 'Reading magnetization components'; call Message("(A)")
  open (UNIT=dataunit,FILE=fname,FORM='UNFORMATTED',status='UNKNOWN')
  read (dataunit) mx
  do i=1,jdim**2
    bx(2*i-1) = mx(i)
    bx(2*i) = 0.0
  end do
  read (dataunit) mx
  do i=1,jdim**2
    by(2*i-1) = mx(i)
    by(2*i) = 0.0
  end do
  read (dataunit) mx
  do i=1,jdim**2
    bz(2*i-1) = mx(i)
    bz(2*i) = 0.0
  end do
  close (unit=dataunit,status='KEEP')

! component-wise FFT of magnetization components
  mess = 'Initializing forward FFT'; call Message("(A)")
  nn(1) = jdim
  nn(2) = jdim
  call fftwnd_f77_create_plan(plan,2,nn,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
  call fftwnd_f77_one(plan,bx,bx)
  call fftwnd_f77_one(plan,by,by)
  call fftwnd_f77_one(plan,bz,bz)
  call fftwnd_f77_destroy_plan(plan)
  mess = 'Magnetization components transformed'; call Message("(A)")

! Compute eqn. 13a-b in the Mansuripur paper.
! The array is numbered in 1-D, although it is actually a 2-D array
! To save space we will use the bx array to store the phase
  d1r=0.0
  d1i=0.0
  d2r=0.0
  d2i=0.0
  d3r=0.0
  d3i=0.0
! loop over y-axis
  do iy=1,jdim
   if (iy.lt.jdim/2) then
    sy=float(iy-1)*fjdim
   else
    sy=float(iy-1)*fjdim-1.0
   endif
! loop over x-axis
   do ix=1,jdim
    iloc=2*(iy-1)*jdim+2*(ix-1)+1
    if (ix.lt.jdim/2) then
     sx=float(ix-1)*fjdim
    else
     sx=float(ix-1)*fjdim-1.0
    endif
! compute normalized frequency components
    s=sqrt(sx**2+sy**2)
    if (s.eq.0.0) then
     sigx=0.0
     sigy=0.0
     izero=iloc
    else 
     sigx=sx/s
     sigy=sy/s
    endif
! compute the products of various vectors
    if (igp) then
      d1r=beam(1)**2*sigx*by(iloc)-beam(2)**2*sigy*bx(iloc)
      d1i=beam(1)**2*sigx*by(iloc+1)-beam(2)**2*sigy*bx(iloc+1)
      d2r=beam(1)*beam(2)*(-bx(iloc)*sigx+by(iloc)*sigy)
      d2i=beam(1)*beam(2)*(-bx(iloc+1)*sigx+by(iloc+1)*sigy)
      d3r=beam(3)*bz(iloc)*(beam(1)*sigy-beam(2)*sigx)
      d3i=beam(3)*bz(iloc+1)*(beam(1)*sigy-beam(2)*sigx)
    end if
    d4r=beam(3)**2*(by(iloc)*sigx-bx(iloc)*sigy)
    d4i=beam(3)**2*(by(iloc+1)*sigx-bx(iloc+1)*sigy)
! compute G_p(ts) (or not)
    gp=1.0
    if (igp) then 
     psig=(beam(1)*sigx+beam(2)*sigy)/beam(3)
     pre=1.0/(psig**2+1.0)/beam(3)**2
     if ((psig.ne.0.0).and.(iloc.ne.izero)) then
      psig=psig*cPi*s*B0t
      gp=pre*sin(psig)/psig
     else
      gp=pre
     end if
    end if
! multiply the whole thing and store in bx
    if (iloc.ne.izero) then 
     bx(iloc)=-gp*(d1i+d2i+d3i+d4i)*B0t/s
     bx(iloc+1)=gp*(d1r+d2r+d3r+d4r)*B0t/s
    else
     bx(iloc)=0.0
     bx(iloc+1)=0.0
    endif
   end do
  end do
  deallocate(by,bz)
! inverse Fourier transform to get the phase
  mess = 'Reciprocal phase computed; starting FFT'; call Message("(A)")
  call fftwnd_f77_create_plan(plan,2,nn,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
  call fftwnd_f77_one(plan,bx,bx)
  call fftwnd_f77_destroy_plan(plan)
  mess = 'FFT complete'; call Message("(A)")

! multiply with prefactor 
  bx = bx*prefac

! save output in file
  mess = 'Store computed phase profile'; call Message("(A)")
  call SafeOpenFile('d1','UNFORMATTED',gname)
  write (dataunit) bx
  call SafeCloseFile('d1','KEEP',gname)
! prepare for TIFF output
  TIFF_nx = jdim
  TIFF_ny = jdim
  TIFF_filename = "phase.tiff"
  do i=1,jdim**2
   mx(i) = bx(2*i-1)
  end do
  mi = minval(mx)
  ma = maxval(mx)
  mx = 255.0*(mx-mi)/(ma-mi)
  write (*,*) 'min(phi), max(phi) = ',mi,ma
! allocate memory for image
  allocate(TIFF_image(0:TIFF_nx-1,0:TIFF_ny-1))

! rotate the calculated array so that it has the proper orientation
  do i=0,TIFF_nx-1
    do j=0,TIFF_ny-1
     iloc = i*TIFF_ny+j+1
     TIFF_image(j,TIFF_ny-1-i) = int(mx(iloc))
    end do
  end do
! create the file
  mess = 'Creating TIFF file: phase.tiff'; call Message("(A)")
  call TIFF_Write_File
  deallocate(TIFF_image)

end subroutine PhaseMap


!
! ###################################################################
! 
!  subroutine CalcLorentzCTF
!
!                                    created: 4/16/97
!  Author: Marc De Graef
!  
!  Description: compute the Lorentz contrast transfer function, which
!               contains defocus, aperture, beam divergence, and astigmatism
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   4/16/97 MDG 1.0 original
!   9/29/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine CalcLorentzCTF(aprad,appos,defocus,defspread,thetac,astig,scl,dimi,dimj)

use local
use constants
use diffraction
use CTF_mod

IMPLICIT NONE

real(kind=sgl)        :: aprad, aprad2, appos(2),defocus,thetac, q, dis,d,scl, &
                         fidim,fjdim,arg,att,t,astig(2),a1,a2,a3,r,plr,defspread
integer(kind=irg)     :: dimi,dimj, ix,iy,i
real(kind=sgl),allocatable      :: idimi(:),jdimj(:)

INTENT(IN)            :: aprad, appos, defocus, thetac, dimi, dimj, astig

  aprad2 = aprad**2
  fidim = 1.0/float(dimi)
  fjdim = 1.0/float(dimj)
! x and y components of spatial frequency vector
  allocate(idimi(dimi),jdimj(dimj))
  idimi = float((/ (i, i=0,dimi-1) /))
  jdimj = float((/ (i, i=0,dimj-1) /))
! wrap around to make (1,1) the origin
  where(idimi.ge.dimi/2) idimi = idimi-float(dimi)
  where(jdimj.ge.dimj/2) jdimj = jdimj-float(dimj)
! apply scaling factors
  idimi = scl*idimi*fidim
  jdimj = scl*jdimj*fjdim
! precalculate constants for transfer function
  a1 = cPi*mLambda
  a2 = (cPi*thetac*defocus)**2
  a3 = 0.5*(cPi*mLambda*defspread)**2
! loop over y axis  
  do iy=1,dimj
! loop over x-axis
    do ix=1,dimi
! spatial frequency
      q=idimi(ix)**2+jdimj(iy)**2
! displaced aperture parameter
      dis=(idimi(ix)-appos(1))**2+(jdimj(iy)-appos(2))**2
! LCTF is allocated in the calling program
      LCTF(ix,iy)=cmplx(0.0,0.0)
! if inside aperture
      if (dis.le.aprad2) then
! polar angle
       r = idimi(ix)/sqrt(q)
       plr = acos(r)
       if ((jdimj(iy).eq.0.0).and.(idimi(ix).lt.0.0)) plr = cPi
       if (jdimj(iy).lt.0.0) plr = 2.0*cPi-plr
       arg=a1*(defocus + astig(1)*cos(2.0*(plr-astig(2))))*q
       att=a2*q + a3*q**2
       if (att.lt.20.0) then
        t=exp(-att)
        LCTF(ix,iy)=cmplx(cos(arg)*t,-sin(arg)*t)
       end if
      end if
    end do
  end do
  deallocate(idimi,jdimj)
end subroutine


end module lorentz_mod

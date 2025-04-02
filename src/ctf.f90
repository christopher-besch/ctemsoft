!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft: ctf.f90                                                             !
! Copyright (c) 1996, 2001, 2002  R.A. Vowels/Marc De Graef/CMU                 !
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
! 
!  FILE: "ctf.f90"
!                                    created: 9/2/01 {11:26:49 AM} 
!                                last update: 9/2/01
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: everything that has to do with contrast transfer
! 
!  History 
!
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   9/2/01 MDG  1.0 f90 version
! ###################################################################
!

module CTF_mod

! the microscope parameters 
type scope_type
  real         :: Cs, &      ! spherical aberration constant [nm]
                  Cc, &      ! chromatic aberration constant [nm]
                  df, &      ! defocus value [nm, under-focus = positive]
                  Ca,phiA, & ! two-fold astigmatism constants [nm, rad]
                  G,phiG,  & ! three-fold astigmatism constants [nm, rad]
                  H,phiH,  & ! four-fold astigmatism constants [nm, rad]
                  K,phiK,  & ! axial coma [nm, rad]
                  dE, &      ! source energy spread [eV]
                  sV, &      ! sigma(V)/V Voltage stability
                  sI, &      ! sigma(I)/I Voltage stability
                  thetac, &  ! beam divergence angle [rad] 
                  V, &       ! accelerating voltage [V]
                  lambda, &  ! electron wavelength [nm]
                  aprad, &   ! diffraction aperture radius [nm^-1]
                  lns, &     ! ln(s), logarithm of signal to noise ratio
                  d(2), &    ! drift vector
                  u(2), &    ! sample vibration vector
                  pixCCD, &  ! pixel size on CCD camera [nm]
                  CCD_N, &   ! number of pixel that are unreliable on CCD (for detector envelope)
                  dq(2), &   ! smallest increment in spatial frequency [nm^-1]
                  qmax, &    ! maximum spatial frequency [nm^-1]
! derived quantities (computed in Set_Scope)
                  sE, &      ! sigma(E)/E Energy stability
                  qzero, &   ! beam divergence angle/wavelength [rad/nm]
                  dfS, &     ! Scherzer defocus [nm]
                  rhoS, &    ! point resolution [nm]
                  dfo, &     ! optimum defocus [nm]
                  delta, &   ! defocus spread [nm]
                  rhoc, &    ! chromatic information limit [nm]
                  rhothetac,&! beam divergence information limit [nm]
                  rhod, &    ! drift information limit [nm]
                  rhou, &    ! sample vibration information limit [nm]
                  rhodet     ! detector information limit [nm]
end type

complex,allocatable   :: ctf(:,:)
real,allocatable      :: qq(:,:),qx(:,:),qy(:,:),polar(:,:)
logical               :: qq_done = .FALSE.

! define a variable of the scope_type
type(scope_type)   :: scope

contains
!
! ###################################################################
! 
!  subroutine Set_Scope 
!
!                                    created: 9/2/01 
!                                last update: 9/2/01
!  Author: Marc De Graef
!  
!  Description: compute all microscope parameters and allocate memory
!               for contrast transfer computation
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   9/ 2/01 MDG 1.0 f90 version
! ###################################################################
subroutine Set_Scope(dimx,dimy,verbose)

use constants
use diffraction
use io

IMPLICIT NONE

integer,intent(IN)    :: dimx,dimy
logical,intent(IN)    :: verbose
real                  :: a

! compute all variables in scope record
 scope%V = mAccVol
 scope%lambda = mLambda
 scope%sE     = scope%dE/scope%V                                        ! variance of energy spread
 scope%delta  = scope%Cc*sqrt(scope%sV**2+4.0*scope%sI**2+scope%sE**2)  ! defocus spread in nm
 scope%qzero  = scope%thetac/scope%lambda                               ! q_0 in nm^-1
 scope%dfS    = sqrt(1.5*scope%Cs*scope%lambda)                         ! Scherzer defocus
 scope%rhoS   = (scope%Cs*scope%lambda**3/6.0)**0.25                    ! Point resolution
 scope%rhoc   = sqrt(cPi*scope%lambda*scope%delta/sqrt(2.0*scope%lns))  ! chromatic information limit
 scope%dfo    = 0.75*scope%Cs*scope%lambda**2/scope%rhoc**2             ! optimum defocus
 a = 0.25
 scope%rhothetac = ((6.0*cPi*a*scope%thetac)/(scope%lambda*sqrt(scope%lns))*scope%rhoS**4)**0.33333 ! beam divergence information limit
 scope%rhod   = cPi*sqrt(scope%d(1)**2+scope%d(2)**2)/sqrt(6.0*scope%lns)             ! drift information limit
 scope%rhou   = cPi*sqrt(scope%u(1)**2+scope%u(2)**2)/sqrt(scope%lns)                 ! vibration information limit
 scope%rhodet = scope%rhoS*((12.0*sqrt(2.0)*cPi*a)/(scope%CCD_N*sqrt(scope%lns)))**0.25 ! detector information limit
! if verbose, then spit it all out 
 if (verbose) then
  mess = 'Defocus spread               [nm] : '; oi_real(1)=scope%delta; call WriteReal(1,'(f8.3)')
  mess = 'Scherzer defocus             [nm] : '; oi_real(1)=scope%dfS  ; call WriteReal(1,'(f8.3)')
  mess = 'Point resolution             [nm] : '; oi_real(1)=scope%rhoS ; call WriteReal(1,'(f8.3)')
  mess = 'Chromatic information limit  [nm] : '; oi_real(1)=scope%rhoc ; call WriteReal(1,'(f8.3)')
  mess = 'Divergence information limit [nm] : '; oi_real(1)=scope%rhothetac ; call WriteReal(1,'(f8.3)')
  mess = 'Drift information limit      [nm] : '; oi_real(1)=scope%rhod ; call WriteReal(1,'(f8.3)')
  mess = 'Vibration information limit  [nm] : '; oi_real(1)=scope%rhou ; call WriteReal(1,'(f8.3)')
  mess = 'Detector information limit   [nm] : '; oi_real(1)=scope%rhodet ; call WriteReal(1,'(f8.3)')
 end if
! allocate memory for contrast transfer function 
 allocate(ctf(0:dimx-1,0:dimy-1),qq(0:dimx-1,0:dimy-1),qx(0:dimx-1,0:dimy-1),qy(0:dimx-1,0:dimy-1),polar(0:dimx-1,0:dimy-1))

end subroutine Set_Scope
!
! ###################################################################
! 
!  subroutine Get_CTF
!
!                                    created: 9/2/01 
!                                last update: 9/2/01
!  Author: Marc De Graef
!  
!  Description: compute contrast transfer function
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   9/ 2/01 MDG 1.0 f90 version
! ###################################################################
subroutine Get_CTF(dimx,dimy)

use constants

IMPLICIT NONE

integer      :: dimx, dimy, i,j
real         :: r,mask(0:dimx-1,0:dimy-1), damp(0:dimx-1,0:dimy-1), &
                arg(0:dimx-1,0:dimy-1), u(0:dimx-1,0:dimy-1)


! reciprocal q-array (only the first time)
 if (.not.qq_done) then 
  do i=0,dimx-1
!  qx(i,0) = float(i-dimx/2)
   qx(i,0) = float(i)-float(i/(dimx/2))*dimx
   do j=1,dimy-1
    qx(i,j) = qx(i,0)
   end do
  end do
  qx = qx*scope%dq(1)
  do j=0,dimy-1
!  qy(0,j) = float(j-dimy/2)
   qy(0,j) = float(j)-float(j/(dimy/2))*dimy
   do i=1,dimx-1
    qy(i,j) = qy(0,j)
   end do
  end do
  qy = qy*scope%dq(2)
  qq = qx**2+qy**2
! polar angle
  do i=0,dimx-1 
   do j=0,dimy-1 
    if (qq(i,j).ne.0.0) then
      r = qx(i,j)/sqrt(qq(i,j))
      polar(i,j) = acos(r)
    end if
   end do
  end do
  where ((qy.eq.0.0).and.(qx.lt.0.0)) polar = cPi
  where(qy.lt.0.0) polar = 2.0*cPi-polar
  polar(dimx/2,dimy/2) = 0.0
  qq_done = .TRUE.
 end if
! aperture function
 mask = 1.0
 where(qq > scope%aprad**2) mask = 0.0
! phase factor
 arg =  cPi*scope%lambda*qq*(scope%df+scope%Ca*cos(2.0*(polar-scope%phiA))) + &
              2.0*cPi*scope%lambda**2*qq**1.5*(scope%G*cos(3.0*(polar-scope%phiG))+3.0*scope%K*cos(polar-scope%phiK))/3.0 + &
              0.5*cPi*scope%lambda**3*qq**2*(scope%H*cos(4.0*(polar-scope%phiH))-scope%Cs)
 ctf = cmplx(cos(arg),sin(arg))
! damping envelope
 u = 1.0/(1.0 + 2.0*(cPi*scope%thetac*scope%delta)**2*qq)
 arg = 0.5*(cPi*scope%lambda*scope%delta*qq)**2 + &
       (cPi*scope%thetac/scope%lambda)**2*(scope%Cs*scope%lambda**3*qq**1.5-scope%df*scope%lambda*sqrt(qq))**2
 ctf = mask * ctf * cmplx(exp(-arg*u))

end subroutine Get_CTF

end module CTF_mod


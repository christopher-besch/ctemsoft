!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:crystal.f90                                                          !
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
! 
!  FILE: "crystal.f90"
!                                    created: 1/5/99 {11:26:49 AM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: everything that has to do with the crystal structure
! 
!  History 
!
!    subroutine CalcMatrices
!    subroutine TransSpace(t,d,inspace,outspace)
!    subroutine TransCoor(t,d,talpha,space)
!    subroutine NormVec(p,space)
!    subroutine CalcCross(p,q,r,inspace,outspace,iv)
!    subroutine MilBrav(p,q,space,d)
!    real function CalcDot(p,q,space)
!    real function CalcLength(p,space)
!    real function CalcAngle(p,q,space)
!
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  1/5/99   MDG 1.0 original
!   7/16/99 MDG 1.1 added error handling and TransCoor
!   4/ 5/00 MDG 1.2 modified TransCoor to include mInvert
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################

module crystal

contains

! ###################################################################
! 
!  subroutine CalcMatrices
!
!  Author: Marc De Graef
!  
!  Description: compute metric tensors and structure matrices
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine CalcMatrices

use local
use error
use constants
use crystalvars

IMPLICIT NONE

real(kind=dbl)     :: det,ca,cb,cg,sa,sb,sg,tg,pirad
integer(kind=irg)  :: i,j


! auxiliary variables for the various tensors
 pirad = cPi/180.0_dbl
 ca = cos(pirad*cell%alpha)
 cb = cos(pirad*cell%beta)
 cg = cos(pirad*cell%gamma)
 sa = sin(pirad*cell%alpha)
 sb = sin(pirad*cell%beta)
 sg = sin(pirad*cell%gamma)
 tg = tan(pirad*cell%gamma)
 if (sg.eq.0.0_dbl) then
! CalcMatrices: invalid gamma angle
  call ErrorMess(3)
  stop
 endif
! define the Kronecker Delta
 cell%krdel = reshape( (/ 1.0_dbl,0.0_dbl,0.0_dbl,0.0_dbl,1.0_dbl,0.0_dbl,0.0_dbl,0.0_dbl,1.0_dbl /), (/3,3/) )
! compute the direct metric tensor [equation 1.5, page 6]
 cell%dmt(1,1) = cell%a**2
 cell%dmt(2,2) = cell%b**2
 cell%dmt(3,3) = cell%c**2
 cell%dmt(1,2) = cell%a*cell%b*cg
 cell%dmt(2,1) = cell%dmt(1,2)
 cell%dmt(1,3) = cell%a*cell%c*cb
 cell%dmt(3,1) = cell%dmt(1,3)
 cell%dmt(2,3) = cell%b*cell%c*ca
 cell%dmt(3,2) = cell%dmt(2,3)
! compute the reciprocal metric tensor as the inverse of the direct
! metric tensor
 det = (cell%a*cell%b*cell%c)**2*(1.0-ca**2-cb**2-cg**2+2.0*ca*cb*cg)
 cell%vol = sqrt(det)
 if (cell%vol.lt.0.1D-6) then
! CalcMatrices: unit cell volume = 0
  call ErrorMess(2)
  stop
 endif
 cell%rmt(1,1) = (cell%b*cell%c*sa)**2
 cell%rmt(2,2) = (cell%a*cell%c*sb)**2
 cell%rmt(3,3) = (cell%a*cell%b*sg)**2
 cell%rmt(1,2) = cell%a*cell%b*cell%c**2*(ca*cb-cg)
 cell%rmt(2,1) = cell%rmt(1,2)
 cell%rmt(1,3) = cell%a*cell%b**2*cell%c*(cg*ca-cb)
 cell%rmt(3,1) = cell%rmt(1,3)
 cell%rmt(2,3) = cell%a**2*cell%b*cell%c*(cb*cg-ca)
 cell%rmt(3,2) = cell%rmt(2,3)
 cell%rmt = cell%rmt/det
! compute the direct structure matrix [equation 1.64, page 57]
 cell%dsm(1,1) = cell%a
 cell%dsm(1,2) = cell%b*cg
 cell%dsm(1,3) = cell%c*cb
 cell%dsm(2,1) = 0.0_dbl
 cell%dsm(2,2) = cell%b*sg
 cell%dsm(2,3) = -cell%c*(cb*cg-ca)/sg
 cell%dsm(3,1) = 0.0_dbl
 cell%dsm(3,2) = 0.0_dbl
 cell%dsm(3,3) = cell%vol/(cell%a*cell%b*sg)
! compute the reciprocal structure matrix [equation 1.65, page 58]
 cell%rsm(1,1) = 1.0_dbl/cell%a
 cell%rsm(1,2) = 0.0_dbl
 cell%rsm(1,3) = 0.0_dbl
 cell%rsm(2,1) = -1.0_dbl/(cell%a*tg)
 cell%rsm(2,2) = 1.0_dbl/(cell%b*sg)
 cell%rsm(2,3) = 0.0_dbl
 cell%rsm(3,1) = cell%b*cell%c*(cg*ca-cb)/(cell%vol*sg)
 cell%rsm(3,2) = cell%a*cell%c*(cb*cg-ca)/(cell%vol*sg)
 cell%rsm(3,3) = (cell%a*cell%b*sg)/cell%vol
end subroutine
!     
! ###################################################################
! 
!  subroutine TransSpace
!
!  Author: Marc De Graef
!  
!  Description: convert vector components from one space (inspace) 
!               to another (outspace); see page 63.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine TransSpace(t,d,inspace,outspace)

use local
use crystalvars
use math

IMPLICIT NONE

real(kind=sgl)  :: t(3), d(3)
character(1)    :: inspace, outspace

intent(IN)      :: t,inspace,outspace
intent(OUT)     :: d

 if (inspace.eq.'d') then
! direct to Cartesian (pre-multiplication)
  if (outspace.eq.'c') then
   d = matmul(cell%dsm,t)
  end if
! direct to reciprocal (post-multiplication)
  if (outspace.eq.'r') then
   d = matmul(t,cell%dmt)
  end if
 end if

 if (inspace.eq.'r') then
! reciprocal to Cartesian (pre-multiplication)
  if (outspace.eq.'c') then
   d = matmul(cell%rsm,t)
  end if
! reciprocal to direct (post-multiplication)
  if (outspace.eq.'d') then
   d = matmul(t,cell%rmt)
  end if
 end if

 if (inspace.eq.'c') then
! Cartesian to direct (post-multiplication)
  if (outspace.eq.'d') then
   d = matmul(cell%rsm,t)
  end if
! Cartesian to reciprocal (post-multiplication)
  if (outspace.eq.'r') then
   d = matmul(t,cell%dsm)
  end if
 end if
end subroutine
!     
! ###################################################################
! 
!  subroutine TransCoor 
!
!  Author: Marc De Graef
!  
!  Description: convert vector components from one reference frame 
!               to another; this is a general coordinate
!               transformation using the old-to-new matrix alpha
!               The details of this routine are summarized in 
!               Table 1.6, page 51, of the textbook. The direction of the 
!               transformation is 'on' (old-to-new) or 'no'
!               (new-to-old).
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  7/16/99 MDG 1.0 original
!  4/ 5/00 MDG 1.1 added support for new mInvert
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine TransCoor(t,d,talpha,space,direction)

use local
use crystalvars
use math

IMPLICIT NONE

real(kind=dbl):: t(3), d(3), talpha(3,3), alinv(3,3)
character(1)  :: space
character(2)  :: direction
logical       :: uni

intent(IN)    :: t,talpha,space,direction
intent(OUT)   :: d

! these matrices are typically unitary, so inverse is simply
! the transpose
 uni = .TRUE.
 if (space.eq.'d') then 
  if (direction.eq.'on') then 
   call mInvert(talpha,alinv,uni)
   d = matmul(t,alinv)
  else
   d = matmul(t,talpha)
  end if
 else
  if (direction.eq.'on') then 
   d = matmul(talpha,t)
  else
   call mInvert(talpha,alinv,uni)
   d = matmul(alinv,t)
  end if
 end if
end subroutine
! ###################################################################
! 
!  function CalcDot
!
!  Author: Marc De Graef
!  
!  Description: computes the dot product between two vectors in
!               real, reciprocal, or Cartesian space; implements
!               equations 1.6 (page 7), and 1.16 (page 15).
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
real function CalcDot(p,q,space)

use local
use crystalvars
use math

IMPLICIT NONE

real(kind=sgl) :: p(3), q(3), x
character(1)   :: space

intent(IN)     :: p,q,space

 if (space.eq.'d') x = dot_product(p,matmul(cell%dmt,q))
 if (space.eq.'r') x = dot_product(p,matmul(cell%rmt,q))
 if (space.eq.'c') x = dot_product(p,q)
 CalcDot = x
end function
! ###################################################################
! 
!  subroutine NormVec  
!
!  Author: Marc De Graef
!  
!  Description: Normalize a vector in real, reciprocal or Cartesian
!               components
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine NormVec(p,space)

use local
use crystalvars

IMPLICIT NONE

real(kind=sgl) :: p(3), x
character(1)   :: space
integer        :: i

intent(IN)     :: space
intent(INOUT)  :: p

 x=CalcLength(p,space)
 if (x.ne.0.0) then 
   p=p/x
 else
   p=(/0.0,0.0,0.0/)
 end if  
end subroutine
! ###################################################################
! 
!  function CalcLength
!
!  Author: Marc De Graef
!  
!  Description: compute the length of a vector in real, reciprocal
!               or Cartesian space
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
real function CalcLength(p,space)

use local
use crystalvars

IMPLICIT NONE

real(kind=sgl) :: p(3), x
character(1)   :: space

intent(IN)     :: p,space

 x = CalcDot(p,p,space)
 CalcLength = sqrt(x)
end function
! ###################################################################
! 
!  function CalcAngle 
!
!  Author: Marc De Graef
!  
!  Description: computes the angle between two vectors in real or
!               reciprocal space;  returned in radians
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
real function CalcAngle(p,q,space)

use local
use crystalvars
use error
use constants

IMPLICIT NONE

real(kind=sgl) :: p(3), q(3), x, y, z, t
character(1)   :: space

intent(IN)     :: p,q,space

 x = CalcDot(p,q,space)
 y = CalcLength(p,space)
 z = CalcLength(q,space)
 if ((y.eq.0.0).or.(z.eq.0.0)) then
! CalcAngle: vector of zero length specified
  call ErrorMess(4)
  stop
 end if
 t = x/(y*z)
 if (t.ge.1.0) then 
  CalcAngle = 0.0
 else 
  if (t.le.-1.0) then 
   CalcAngle = cPi
  else 
   CalcAngle = acos(t)
  end if
 end if
end function
! ###################################################################
! 
!  subroutine CalcCross
!
!  Author: Marc De Graef
!  
!  Description: computes the cross product between two vectors and
!               expresses it in either real space or reciprocal space.
!               The output can also be expressed in the standard
!               Cartesian reference frame.  The switch iv indicates
!               whether the result should be scaled by the unit cell
!               volume. More information in section 1.3.5, page 18.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine CalcCross(p,q,r,inspace,outspace,iv)

use local
use crystalvars
use math

IMPLICIT NONE 

real(kind=sgl)    :: p(3), q(3), r(3), x(3), vl
character(1)      :: inspace, outspace
integer(kind=irg) :: iv

intent(IN)        :: p,q,inspace,outspace,iv
intent(OUT)       :: r

 if (iv.eq.1) then 
  vl = cell%vol
 else
  vl = 1.0
 endif
 if (inspace.eq.'d') then
  r(1) = vl*(p(2)*q(3)-p(3)*q(2))
  r(2) = vl*(p(3)*q(1)-p(1)*q(3))
  r(3) = vl*(p(1)*q(2)-p(2)*q(1))
  if (outspace.eq.'d') then 
   x = matmul(r,cell%rmt)
   r = x
  end if
  if (outspace.eq.'c') then 
   x = matmul(cell%rsm,r)
   r = x
  end if
 end if
 if (inspace.eq.'r') then
  r(1) = (p(2)*q(3)-p(3)*q(2))/vl
  r(2) = (p(3)*q(1)-p(1)*q(3))/vl
  r(3) = (p(1)*q(2)-p(2)*q(1))/vl
  if (outspace.eq.'r') then 
   x = matmul(r,cell%dmt)
   r = x
  end if
  if (outspace.eq.'c') then 
   x = matmul(cell%dsm,r)
   r = x
  end if
 end if
 if (inspace.eq.'c') then
  r(1) = p(2)*q(3)-p(3)*q(2)
  r(2) = p(3)*q(1)-p(1)*q(3)
  r(3) = p(1)*q(2)-p(2)*q(1)
 end if
end subroutine
! ###################################################################
! 
!  subroutine MilBrav
!
!  Author: Marc De Graef
!  
!  Description: conversion from Miller to Miller-Bravais indices for
!               directions.  The switch d is either '34' or '43'.
!               implements equations 1.31 and 1.32, pages 24-25.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  11/15/00 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine MilBrav(p,q,d)

use local

IMPLICIT NONE

integer(kind=irg)      :: p(3), q(4), i, j
real(kind=sgl)         :: r(4), rm
character(2)           :: d

intent(IN)             :: d
intent(INOUT)          :: p,q

 if (d.eq.'43') then 
! equation 1.31
! these will always be integers, so no reduction is required
  p(1) = q(1)-q(3)
  p(2) = q(2)-q(3)
  p(3) = q(4)
 else
! equation 1.32
! there is no need to divide by 3, since that would be taken out 
! by the reduction to integers in the next step
  r(1) = float(2*p(1)-p(2))
  r(2) = float(2*p(2)-p(1))
  r(3) = -float(p(1)+p(2))
  r(4) = float(3*p(3))
! next reduce to common integers
! first, find the non-zero minimum index
  rm = 100.0
  do i=1,4 
   if ((abs(r(i)).lt.rm).and.(r(i).gt.0.0)) then
    rm = abs(r(i))
   end if
  end do
! then check if this index is a common divider of the others
  j = 0
  do i=1,4
   r(i) = r(i)/rm
   if ((r(i)-mod(r(i),1.0)).eq.0.0) j=j+1
  end do
  if (j.eq.4) rm=1.0
   q = int(r*rm)
 end if
end subroutine
! ###################################################################
! 
!  subroutine GetLatParm
!
!                                    created: 10/13/98 {9:29:46 AM} 
! 
!  Author: Marc De Graef
!  
!  Description: asks for crystal system and lattice parameters
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine GetLatParm

use local
use io
use symmetryvars
use crystalvars

IMPLICIT NONE

integer(kind=irg)        :: i

 mess = ' Select the crystal system : '; call Message("(A)")
 mess = '  1. Cubic '; call Message("(A)")
 mess = '  2. Tetragonal '; call Message("(A)")
 mess = '  3. Orthorhombic '; call Message("(A)")
 mess = '  4. Hexagonal '; call Message("(A)")
 mess = '  5. Trigonal '; call Message("(A)")
 mess = '  6. Monoclinic '; call Message("(A)")
 mess = '  7. Triclinic '; call Message("(A)")
 mess = ' '; call Message("(A)")
 mess = ' Note about the trigonal system:'; call Message("(A)")
 mess = ' -------------------------------'; call Message("(A)")
 mess = ' Primitive trigonal crystals are defined'; call Message("(A)")
 mess = ' with respect to a HEXAGONAL reference frame.'; call Message("(A)")
 mess = ' Rhombohedral crystals can be referenced with'; call Message("(A)")
 mess = ' respect to a HEXAGONAL basis (first setting),'; call Message("(A)")
 mess = ' or with respect to a RHOMBOHEDRAL basis (second'; call Message("(A)")
 mess = ' setting).  The default setting for trigonal'; call Message("(A)")
 mess = ' symmetry is the hexagonal setting.  When you '; call Message("(A)")
 mess = ' select crystal system 5 above, you will be '; call Message("(A)")
 mess = ' prompted for the setting.'; call Message("(A)")
 mess = ' '; call Message("(A)")
 mess = ' crystal system ---> '; call GetInt(1)
 cell%xtal_system = io_int(1)
! make sure the symmetry operations will be reduced to the 
! fundamental unit cell
 cell%SYM_reduce=.TRUE.
 hexset=.FALSE.
! deal with the rhombohedral vs. hexagonal setting
! (the rhombohedral axes are considered as the second setting)
 i=-1
 cell%SYM_trigonal=.FALSE.
 cell%SYM_second=.FALSE.
 if (cell%xtal_system.eq.5) then
  cell%SYM_trigonal=.TRUE.
  mess = 'Enter 1 for rhombohedral lattice parameters,'; call Message("(A)")
  mess = '0 for hexagonal lattice parameters :'; call GetInt(1)
  i = io_int(1)
  if (i.eq.0) then
   cell%xtal_system=4
  else
   cell%SYM_second=.TRUE.
  end if
 end if
 mess = 'Enter lattice parameters'; call Message("(A)")
! put default values based on cubic symmetry, then change them later
 mess = '    a [nm] = '; call GetReal(1)
 cell%a = io_real(1)
 cell%b=cell%a 
 cell%c=cell%a 
 cell%alpha=90.0
 cell%beta=90.0
 cell%gamma=90.0
! now get the proper lattice parameters
 select case (cell%xtal_system)
  case (1)
! tetragonal
  case (2)
   mess = '    c [nm] = '; call GetReal(1)
   cell%c = io_real(1)
! orthorhombic
  case (3)
   mess = '    b [nm] = '; call GetReal(1)
   cell%b = io_real(1)
   mess = '    c [nm] = '; call GetReal(1)
   cell %c = io_real(1)
! hexagonal
  case (4)
   mess = '    c [nm] = '; call GetReal(1)
   cell%c = io_real(1)
   cell%gamma=120.0
! rhombohedral 
  case (5)
   mess = '    alpha [deg] = '; call GetReal(1)
   cell%alpha = io_real(1)
   cell%beta=cell%alpha
   cell%gamma=cell%alpha
! monoclinic   
  case (6)
   mess = '    b [nm] = '; call GetReal(1)
   cell%b = io_real(1)
   mess = '    c [nm] = '; call GetReal(1)
   cell%c = io_real(1)
   mess = '    beta  [deg] = '; call GetReal(1)
   cell%beta = io_real(1)
! triclinic    
  case (7) 
   mess = '    b [nm] = '; call GetReal(1)
   cell%b = io_real(1)
   mess = '    c [nm] = '; call GetReal(1)
   cell%c = io_real(1)
   mess = '    alpha [deg] = '; call GetReal(1)
   cell%alpha = io_real(1)
   mess = '    beta  [deg] = '; call GetReal(1)
   cell%beta = io_real(1)
   mess = '    gamma [deg] = '; call GetReal(1)
   cell%gamma = io_real(1)
 end select
! if trigonal symmetry was selected in the first setting,
! then the xtal_system must be reset to 5
 if (cell%SYM_trigonal) then
  cell%xtal_system=5
 end if
! if hexagonal setting is used, then Miller-Bravais indices must be enabled
 if ((cell%xtal_system.eq.4).OR.((cell%xtal_system.eq.5).AND.(.not.cell%SYM_second))) then
  hexset = .TRUE.
 else 
  hexset = .FALSE.
 end if
end subroutine
!
! ###################################################################
! 
!  subroutine GetAsymPos
!
!                                    created: 10/13/98 {9:29:46 AM} 
! 
!  Author: Marc De Graef
!  
!  Description: get the fractional coordinates for all atoms
!               in the asymmetric unit
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine GetAsymPos

use local
use io
use crystalvars

IMPLICIT NONE

logical                 :: more
character(1)            :: ans,list(256)
real(kind=sgl)          :: pt(5)
integer(kind=irg)       :: j

 more=.TRUE.
 cell%ATOM_ntype = 0
 mess = ' Enter atoms in asymmetric unit '; call Message("(A)")
 call DisplayElements
 do while (more)
  cell%ATOM_ntype = cell%ATOM_ntype + 1
  write (*,"(' ->  Atomic number : ',$)")
  read (*,"(I2)") cell%ATOM_type(cell%ATOM_ntype)
  write (*,"(' ->  Fractional coordinates, site occupation, and Debye-Waller Factor [nm^2] : ')") 
  do j=1,256
   list(j)(1:1) = ' '
  end do
  read (*,fmt="(256A)") list
  call extractposition(list,pt) 
  cell%ATOM_pos(cell%ATOM_ntype,1:5) = pt(1:5)
  write (*,"(1x,4(F10.7,2x),F10.7)")  (cell%ATOM_pos(cell%ATOM_ntype,j),j=1,5)
  write (*,"(' ->  Another atom ? (y/n) ',$)")
  read (*,"(A1)") ans
  if ((ans.eq.'y').or.(ans.eq.'Y')) then 
   more=.TRUE.
  else
   more=.FALSE.
  end if 
 end do
end subroutine
! ###################################################################
! 
!  subroutine DisplayElements
!
!                                    created: 10/13/98 {9:29:46 AM}
!
!  Author: Marc De Graef
!  
!  Description: Display the elements of the periodic table
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine DisplayElements

use local
use io
use constants

IMPLICIT NONE

integer(kind=irg) ::  i

 mess = '-------Periodic Table of the Elements--------'; call Message("(/A/)")
 do i=1,98
  if (mod(i,7).eq.0) then
   write (*,"(1x,i3,':',A2)") i,ATOM_sym(i)
  else
   write (*,"(1x,i3,':',A2,2x,$)") i,ATOM_sym(i)
  end if  
 end do
 mess = '----------------------------------------------'; call Message("(A/)")
end subroutine
! ###################################################################
! 
!  subroutine extractposition
!
!                                    created: 10/13/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: take a string and extract the atom positions, Debye-Waller
!               factor and occupation parameter.  The atom positions may
!               be entered as real numbers, or as fractions.  There are probably 
!               other (and maybe better or shorter) ways of doing this, but
!               since it works, why change it ... ?
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   5/25/01 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine extractposition(list,pt)

use local

IMPLICIT NONE

character(1)                :: list(256)
integer(kind=irg)           :: comma(6),slash(5),period(5),ccnt,scnt,pcnt,pp,i,j,hcnt, &
                               ip,ipt,icnt,nd,n,k,ns
integer(kind=irg),parameter :: nmb(48:57)=(/0,1,2,3,4,5,6,7,8,9/)
real(kind=dbl)              :: nominator,denominator,x
real(kind=sgl)              :: pt(5)
logical                     :: hasperiod

 comma(1:6) = 0
 slash(1:5) = 0
 period(1:5) = 0
 ccnt = 0
 scnt = 0
 pcnt = 0
 j = 0
 hcnt = 0
! count characters and search for , . and /
 ccnt = ccnt+1
 comma(ccnt) = 0
 do i=1,256
  if (list(i)(1:1).ne.' ') j=j+1
  if (list(i)(1:1).eq.',') then 
   ccnt = ccnt+1
   comma(ccnt)=i
  end if
  if (list(i)(1:1).eq.'/') then 
   scnt = scnt+1
   slash(scnt)=i
  end if
  if (list(i)(1:1).eq.'.') then 
   pcnt = pcnt+1
   period(pcnt)=i
  end if
 end do 
 ccnt = ccnt+1
 comma(ccnt) = j+1
 do while (ccnt.lt.6) 
  ccnt = ccnt+1
  comma(ccnt) = comma(ccnt-1)+1
 end do
! interpret the string
 j = 1
 ip = 1
 icnt = 0
 ipt = 1
 pp = 1
 do i=1,ccnt-1
! is it a real number or a fraction ?
  if (((slash(j).lt.comma(i+1)).and.(scnt.gt.0)).and.(j.le.scnt)) then
! it is a fraction;  get the nominator
   nd = slash(j)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   nominator = dble(n)
   ip = slash(j)+1
! and then the denominator
   nd = comma(i+1)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   denominator = dble(n)
! and fill in the entire range
   pt(ipt) = sngl(nominator/denominator)
   ipt = ipt+1
   ip = comma(i+1)+1
   j=j+1
  else
! no, it is a real number, possibly without a period
! is there a period in this number ?
   if ((period(pp).gt.comma(i)).and.(period(pp).lt.comma(i+1))) then
     hasperiod = .TRUE.
   else
     hasperiod = .FALSE.
   endif
   nd = comma(i+1)-ip
   if (hasperiod) then 
    if (period(pp).eq.comma(i)+1) then
     x = 0.D0
     ns = 2
    else
     x = dble(nmb(ichar(list(ip)(1:1))))
     ns = 3
    end if 
    do k=ns,nd
     x = x + 10.D0**(ns-k-1)*dble(nmb(ichar(list(ip+k-1)(1:1))))
    end do
    pt(ipt)= sngl(x)
    ipt=ipt+1
    ip = comma(i+1)+1
    pp = pp+1
   else
    nd = comma(i+1)-ip
    n = 0
    do k=0,nd-1
     n = 10*n+nmb(ichar(list(ip+k)(1:1)))
    end do
    pt(ipt) = float(n)
    ipt=ipt+1
    ip = comma(i+1)+1
   end if
  end if
 end do 
! set default values
 if (pt(4).eq.0.0) pt(4) = 1.0
end subroutine
! ###################################################################
!
!  subroutine GetOR
!
!                                    created: 10/20/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: ask for orientation relation between two crystals
! 
!  History 
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine GetOR(orel)

use local
use crystalvars
use io

IMPLICIT NONE

type(orientation) :: orel
real(kind=sgl)    :: c1,c2
integer(kind=irg) :: i

intent(OUT)       :: orel

 c1 = 1.0
 c2 = 1.0
 do while ((c1.ne.0.0).or.(c2.ne.0.0))
  mess = 'Enter orientation relation in following form:'; call Message("(A)")
  mess = 'planes:     h_A,k_A,l_A,h_B,k_B,l_B '; call Message("(A)")
  mess = 'directions: u_A,v_A,w_A,u_B,v_B,w_B '; call Message("(A)")
  mess = 'Plane normals :'; call GetInt(6); 
  orel%gA(1:3) = float(io_int(1:3))
  orel%gB(1:3) = float(io_int(4:6))
  mess = 'Directions    :'; call GetInt(6); 
  orel%tA(1:3) = float(io_int(1:3))
  orel%tB(1:3) = float(io_int(4:6))
! check for orthonormality using zone equation
  c1=sum(orel%tA*orel%gA)
  if (c1.ne.0.0) then
   mess = 'Plane does not contain direction (crystal A)'; call Message("(A)")
  end if
  c2=sum(orel%tB*orel%gB)
  if (c2.ne.0.0) then
   mess = 'Plane does not contain direction (crystal B)'; call Message("(A)")
  end if
 end do
end subroutine
end module

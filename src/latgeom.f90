!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft: latgeom.f90                                                         !
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
!  None - none
! 
!  FILE: "latgeom.f90"
!                                    created: 1/13/99 {11:26:49 AM} 
!                                 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: general lattice geometry computations
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  1/5/99   MDG 1.0 original
!  5/21/01  MDG 2.0 f90
! ###################################################################
!
program latgeom

use local
use io
use files
use crystalvars
use crystal
use constants

integer(kind=irg)        :: isel,another,iv(3)
real(kind=sgl)           :: v1(3),v2(3),vc(3),phi,p,q,r
character(1)             :: sp,sp2

 progname = 'latgeom.f90'
 progdesc = 'Simple lattice geometry program'
 call CTEMsoft
 
! load crystal structure data
 call CrystalData
! set up loop
 another=1
 do while (another.eq.1) 
  call PrintMenu(isel)
  sp='d'
  if (mod(isel,2).eq.0) sp='r'
  if (isel.lt.3) then
   mess = 'Enter vector components : '; call GetInt(3);
   do i=1,3
    v1(i) = float(io_int(i))
   end do
   p=CalcLength(v1,sp)
   if (sp.eq.'d') then 
    mess = '-> Length [nm] = '
   else
    mess = '-> Length [nm-1] = '
   end if
   oi_real(1) = p; call WriteReal(1,"(2x,F10.6)")
  else
   mess = 'Enter first vector components  : '; call GetInt(3);
   do i=1,3
    v1(i) = float(io_int(i))
   end do
   mess = 'Enter second vector components : '; call GetInt(3);
   do i=1,3
    v2(i) = float(io_int(i))
   end do
   if (isel.lt.5) then 
    p=CalcLength(v1,sp)
    q=CalcLength(v2,sp)
    r=CalcDot(v1,v2,sp)
    phi=CalcAngle(v1,v2,sp)*180.0/cPi
    mess = '-> Angle [deg] = '; oi_real(1)=phi; call WriteReal(1,"(2x,F8.4)")
   else
    if (sp.eq.'d') sp2='r'
    if (sp.eq.'r') sp2='d'
    call CalcCross(v1,v2,vc,sp,sp2,0)
    do i=1,3 
       oi_int(i)=int(vc(i))
    end do
    if (sp.eq.'d') then
     call WriteInt(3,"('p x q = (',2(i3,1x),i3,')')")
    else
     call WriteInt(3,"('p x q = [',2(i3,1x),i3,']')")
    end if
   end if
  end if
  mess = 'Another computation? (1/0) : '; call GetInt(1); another = io_int(1)
 end do 
end program
!
!
!
subroutine PrintMenu(isel)

use local
use io
        
integer  :: isel

 mess = 'Select from the following options: '; call Message("(A)")
 mess = ' '; call Message("(A)")
 mess = '[1] length of direct space vector'; call Message("(A)")
 mess = '[2] length of reciprocal space vector'; call Message("(A)")
 mess = '[3] angle between direct space vectors'; call Message("(A)")
 mess = '[4] angle between reciprocal space vectors'; call Message("(A)")
 mess = '[5] cross product, direct space vectors'; call Message("(A)")
 mess = '[6] cross product, reciprocal space vectors'; call Message("(A)")
 mess = ' '; call Message("(A)")
 mess = 'Enter selection: '; call GetInt(1); isel = io_int(1)
end subroutine

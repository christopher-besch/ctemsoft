!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft: orbit.f90                                                           !
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
!  FILE: "orbit.f90"
!                                    created: 1/5/99 {11:26:55 AM}
!                                
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: print orbit of an atom position for a given space group
!               (all atom positions are reduced to the basic unit cell)
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  1/5/99   MDG 1.0 original
!  5/20/01  MDG 2.0 f90
! ###################################################################
program orbit

use local
use crystalvars
use symmetryvars
use symmetry
use files
use io

real(kind=dbl)           :: ctmp(192,3)
integer(kind=irg)        :: i,m,n,ans

 progname = 'orbit.f90'
 progdesc = 'List the orbit of a given position'
 call CTEMsoft
 
 SG % SYM_reduce=.TRUE.
 call CrystalData
 ans = 1
 do while (ans.eq.1)
  mess = ' Enter fractional coordinates of atom : '; call GetReal(3)
  cell % ATOM_pos(1,1:3) = io_real(1:3)
  m=1
  call CalcOrbit(m,n,ctmp)
  mess = ' # equivalent atoms in orbit          : '; oi_int(1) = n; call WriteInt(1,"(I3)")
! spit them out 
  do i=1,n
   write (*,"(1x,i3,' -> (',f7.4,',',f7.4,',',f7.4,');   ',$)") i,(ctmp(i,j),j=1,3)
   if (mod(i+1,2).eq.1) write (*,"(A1)") ' '
  end do
  mess = 'Another orbit ? (1/0) '; call GetInt(1)
  ans = io_int(1)
 end do 

end program

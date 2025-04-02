!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft: listSG.f90                                                         !
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
!  FILE: "listSG.f90"
!                                    created: 1/5/99 {11:26:49 AM}
!                                last update: 6/5/2001 {7:33:35 PM}
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: list general equivalent positions for a 
!               given space group.  This is a simple program
!               to illustrate how one can use the space group matrices
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  1/5/99   MDG 1.0 original
!  5/29/01  MDG 2.0 f90
! ###################################################################
program listSG

use local
use io
use crystalvars
use symmetryvars
use symmetry

IMPLICIT NONE

character(3)   :: pos
character(1)   :: sgn
integer        :: p(4),ii,jj,i
real           :: ppp

 progname = 'listSG.f90'
 progdesc = 'List equivalent positions for arbitrary space group'
 call CTEMsoft

cell % SYM_reduce=.TRUE.
pos = 'xyz'
mess = 'Enter Space Group number : '; 
call GetInt(1); 
cell % SYM_SGnum = io_int(1)
call GenerateSymmetry(.TRUE.)

! loop over all symmetry matrices
 mess = 'Space Group Symbol       : '
 call WriteStr(SYM_SGname(cell % SYM_SGnum))
 mess = 'number of operators      : '; oi_int(1) = SG % SYM_MATnum
 call WriteInt(1,"(I3)")
 do i=1,SG % SYM_MATnum
  oi_int(1) = i; 
  call WriteInt(1,"(1x,i3,2x,'-> (',$)")
! loop over all rows
  do ii=1,3
! loop over colums (get the numbers)
   do jj=1,3 
    p(jj)=int(SG % SYM_data(i,ii,jj)) 
   end do
   ppp = sngl(SG % SYM_data(i,ii,4))
! print each entry 
! first the part containing x, y, and/or z
   do jj=1,3
    if (p(jj).ne.0) then
     mess(1:1)='+'
     if (p(jj).eq.-1) mess(1:1)='-'
     mess(2:2) = pos(jj:jj)
     call Message("(A2,$)")
    end if
   end do 
! if there is a translation component, print it
   if (ppp.ne.0.0) then
    mess(1:1)='+'
    if (ppp.lt.0.0) mess(1:1)='-'
    call Message("(A1,$)");
    oi_real(1) = abs(ppp); 
    call WriteReal(1,"(f5.3,$)")
   end if
! print a comma, or close the brackets and do a newline
   if (ii.ne.3) then 
     mess(1:1) = ','; call Message("(A1,$)")
    else
     mess(1:1) = ')'; call Message("(A1)")
   end if
  end do
 end do 
end program

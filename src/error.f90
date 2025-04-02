!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:error.f90                                                            !
! Copyright (c) 2001, 2002  Marc De Graef/Carnegie Mellon University            !
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
!  FILE: "error.f90"
!                                    created: 1/5/99 {11:26:49 AM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: everything that has to do with error handling
! 
!  History 
!
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   1/5/99  MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################

module error

contains

!
! ###################################################################
!
!  subroutine FatalError
!
!  Author: Marc De Graef
! 
!  Description: write an error message to the screen
!               and exit from the current program
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine FatalError(var1,var2)

use local
use io

IMPLICIT NONE

character(*)  :: var1,var2

intent(IN)    :: var1,var2

 mess = mess//var1//var2; call Message("(' Fatal error in routine ',A)"); stop
end subroutine

!
! ###################################################################
!
!  subroutine ErrorMess
!
!  Author: Marc De Graef
! 
!  Description: write an error message to the screen
!               and exit from the current program
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################

subroutine ErrorMess(i)

use local
use io

IMPLICIT NONE

integer(kind=irg)       :: i
character(60),parameter :: errors(10) = &
     (/"mInvert: transformation matrix zero determinant             ", &  !  1
       "CalcMatrices: unit cell volume = 0                          ", &  !  2
       "CalcMatrices: invalid gamma angle                           ", &  !  3
       "CalcAngle: vector of zero length specified                  ", &  !  4
       "ShortestG: beam direction cannot be [0,0,0]                 ", &  !  5
       "RankReflections: variable numr should be increased          ", &  !  6
       "DumpZAP: maximum number of reflections per pattern exceeded ", &  !  7
       "BWshow: unkown data format                                  ", &  !  8
       "                                                            ", &  !  9
       "                                                            "/)
intent(IN)              :: i

 mess = errors(i); call Message("(' Fatal error in routine ',A)"); stop
end subroutine

end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:local.f90                                                            !
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
!  FILE: "local.f90"
!                                    created: 1/8/99 {4:38:45 PM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: local definitions of single and double precision,
!               general constants and variables
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  1/8/99   MDG 1.0 original
!  5/6/01   MDG 2.0 f90
! 11/27/01  MDG 2.1 added sgl and dbl kinds
! 12/08/01  MDG 2.2 added CTEMsoft subroutine
! ###################################################################

module local

! this module must be "use"d by every program, subroutine, and function !

! The entire CTEMsoft package should be processor independent.  This can
! be accomplished by the use of the "kind" parameters.
! Define the "kind" parameters for single and double precision reals,
  integer,parameter        :: sgl = SELECTED_REAL_KIND(p=6,r=37)
  integer,parameter        :: dbl = SELECTED_REAL_KIND(p=13,r=200)
! Define the "kind" parameters for short and regular integers,
  integer,parameter        :: ish = SELECTED_INT_KIND(3)  ! integer short
  integer,parameter        :: irg = SELECTED_INT_KIND(9)  ! integer regular
! source code version number
  character(8), parameter  :: scversion="2.0/2002"
! author information
  character(13), parameter :: username="Marc De Graef"
  character(26), parameter :: userlocn="Carnegie Mellon University"
! program descriptor (to be defined by each individual program)
  character(15)            :: progname="default"
  character(60)            :: progdesc="default"
! input/output information
! psunit    = Postscript output unit number
! dataunit  = Data unit number (for *.xtal files and other in/output)
  integer, parameter       :: psunit=20, dataunit=21, dataunit2=22, dataunit3=23
! parameters governing the size of varous arrays
! maxpasym  = Maximum number of positions in asymmetric unit
  integer, parameter       :: maxpasym=50
! strucdef is a logical variable to determine whether or not a structure has been loaded;
! hexset determines whether or not 4-index hexagonal indices should be used.
  logical                  :: strucdef,hexset
! where is the postscript viewer on this system ?
  character(18),parameter  :: psviewer="/usr/local/bin/gv "   ! note the space at the end !!!
contains

! ###################################################################
! 
!  subroutine CTEMsoft
!
!  Author: Marc De Graef
!  
!  Description: Prints some general information about the program
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  12/08/01 MDG 1.0 original
! ###################################################################
subroutine CTEMsoft

 write (*,"(//1x,'CTEMsoft version ',A8,', Copyright (C) 2001, 2002 Marc De Graef/CMU')") scversion
 write (*,"(1x,'CTEMsoft comes with ABSOLUTELY NO WARRANTY; for details type `GPL''.')")
 write (*,"(1x,'This is free software, and you are welcome to redistribute it')")
 write (*,"(1x,'under certain conditions; type `GPL'' for details.'//)")

 write (*,"(1x,'Program name : ',A15)") progname
 write (*,"(1x,A60//)") progdesc

end subroutine CTEMsoft

end module

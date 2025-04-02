!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:io.f90                                                               !
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
! -*-Fortran-*-
! ###################################################################
! 
!  FILE: "io.f90"
!                                    created: 1/5/99 {11:26:49 AM} 
!                                last update: 6/1/2001 {6:51:54 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: everything that has to do with input-output
! 
!  History 
!
! 
!  modified by  rev reason
!  -------- --- --- -----------
!    1/5/99 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
!

module io

! general purpose variables for all IO

integer,allocatable         :: io_int(:)
real,allocatable            :: io_real(:)
integer                     :: oi_int(100)
real                        :: oi_real(100)
character(80)               :: mess

contains
!
! ###################################################################
! 
!  subroutine GetStr    
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read a string from standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
subroutine GetStr(oi_str,dim)

integer          :: dim
character(dim)   :: oi_str

 oi_str(1:dim) = " "
 write (*,fmt="(' ',A,' ',$)") trim(mess)
 read (*,fmt="(A)") oi_str
end subroutine
!
! ###################################################################
! 
!  subroutine GetInt
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read one or more integers from input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
subroutine GetInt(dim)

integer     :: i,dim

! first, assume that there is nothing usefull in io_int and deallocate it
 if (allocated(io_int)) deallocate(io_int)
! then initialize it with the proper dimensions
 allocate(io_int(dim))
 write (*,fmt="(' ',A,' ',$)") trim(mess)
 read (*,fmt="(100I8)") (io_int(i),i=1,dim)
 mess(1:80) = ' '
end subroutine
!
! ###################################################################
! 
!  subroutine GetReal
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read one or more reals from input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
subroutine GetReal(dim)

integer     :: dim

! first, assume that there is nothing usefull in io_real and deallocate it
 if (allocated(io_real)) deallocate(io_real)
! then initialize it with the proper dimensions
 allocate(io_real(dim))
 write (*,fmt="(' ',A,' ',$)") trim(mess)
 read (*,fmt="(100F16.8)") (io_real(i),i=1,dim)
 mess(1:80) = ' '
end subroutine
! ###################################################################
! 
! function Message
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 12/14/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: dump a message on standard output
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
subroutine Message(frm)

character(*)  :: frm

 write (*,fmt=frm) trim(mess)
 mess(1:80) = ' '
end subroutine
! ###################################################################
! 
! function WriteInt
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 12/14/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: write integer data to standard output
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
subroutine WriteInt(dim,frm)

character(*)    :: frm
integer         :: i,dim

 write (*,fmt="(' ',A,$)") trim(mess)
 write (*,fmt=frm) (oi_int(i),i=1,dim)
 mess(1:80) = ' '
end subroutine
! ###################################################################
! 
! function WriteReal
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 12/14/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: write real data to standard output
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
subroutine WriteReal(dim,frm)

character(*)    :: frm
integer         :: i,dim

 write (*,fmt="(' ',A,$)") trim(mess)
 write (*,fmt=frm) (oi_real(i),i=1,dim)
 mess(1:80) = ' '
end subroutine
! ###################################################################
! 
! function WriteStr
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 12/14/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: write string data to standard output
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
subroutine WriteStr(m)

character(*)   ::  m

 write (*,fmt="(' ',A,A)") trim(mess),m
end subroutine
! ###################################################################
! 
!  subroutine GetCameraLength
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 10/20/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: ask for camera length, in mm
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
subroutine GetCameraLength(camlen)

real :: camlen

 mess = 'Enter the diffraction camera length L [mm]: '; call GetInt(1)
 camlen = float(io_int(1))
end subroutine


end module

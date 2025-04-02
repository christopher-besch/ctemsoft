!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:  timing.f90                                                         !
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
! -Fortran-90 Source Code-
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "timing.f90"
!                                    created: 4/28/01  {9:29:46 AM} 
!                                 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: timing routines
!               
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
! 11/19/01  MDG 1.0 f90
! ###################################################################
! 
! 

module timing

use local
use io

IMPLICIT NONE

real(kind=sgl)      :: TIME_t_count,TIME_unit_count,TIME_interval,TIME_fraction
integer(kind=irg)   :: TIME_newcount,TIME_count_rate,TIME_count_max,TIME_count


contains

! ###################################################################
! 
!  subroutine Time_reset
!
!  Author: Marc De Graef
!  
!  Description: reset time recording
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
! ###################################################################
subroutine Time_reset

TIME_t_count = 0.0
TIME_unit_count = 0.0
TIME_count = 0
TIME_newcount = 0
TIME_count_rate = 0
TIME_count_max = HUGE(0)

end subroutine Time_reset
! ###################################################################
! 
!  subroutine Time_report
!
!  Author: Marc De Graef
!  
!  Description: how often should the timing information be shown ?
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
! ###################################################################
subroutine Time_report(interval)

use local
use io

IMPLICIT NONE

integer(kind=irg),intent(IN)   :: interval

  TIME_interval = 0.01*float(interval)
  TIME_fraction = TIME_interval
  mess = 'Starting computation'; call Message("(/A)")

end subroutine Time_report
! ###################################################################
! 
!  subroutine Time_start
!
!  Author: Marc De Graef
!  
!  Description: start time recording
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
! ###################################################################
subroutine Time_start

! start the timing of the computation
 call Time_reset
 call system_clock(TIME_count,TIME_count_rate,TIME_count_max)

end subroutine Time_start
! ###################################################################
! 
!  subroutine Time_estimate
!
!  Author: Marc De Graef
!  
!  Description: estimate duration of a single computation step
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
! ###################################################################
subroutine Time_estimate(numk)

use local
use io

IMPLICIT NONE

integer,intent(IN)     :: numk
integer(kind=irg)      :: TIME_nc

! get the current time
 call system_clock(TIME_nc,TIME_count_rate,TIME_count_max)
 TIME_newcount = TIME_nc
 TIME_t_count = float(TIME_newcount-TIME_count)/float(TIME_count_rate)
 TIME_unit_count = TIME_t_count
 mess = ' Time for first computation step [s, typically overestimate] :';
 oi_real(1) = TIME_unit_count
 call WriteReal(1,"(F10.5)")
 mess = '  Anticipated total computation time :'; call Message("(A$)")
 call PrintTime(TIME_unit_count*float(numk))
! 
end subroutine Time_estimate
! ###################################################################
! 
!  subroutine Time_remaining
!
!  Author: Marc De Graef
!  
!  Description: estimate the remaining time
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
! ###################################################################
subroutine Time_remaining(ik,numk)

use local
use io

IMPLICIT NONE

integer,intent(IN)   :: ik,numk
integer(kind=irg)    :: TIME_nc


 mess = ' '; oi_int(1) = nint(100.0*TIME_fraction); call WriteInt (1,"(1x,I2,' % completed')") 
 TIME_fraction = TIME_fraction + TIME_interval
! get the current time
 call system_clock(TIME_nc,TIME_count_rate,TIME_count_max)
 TIME_newcount = TIME_nc
! and print it
 mess = ' Total computation time [s] ' 
 TIME_t_count = float(TIME_newcount-TIME_count)/float(TIME_count_rate)
 oi_real(1) = TIME_t_count
 call WriteReal(1,"(F)")
! reset the time per unit
 TIME_unit_count = TIME_t_count/float(ik)
! print estimated remaining time
 mess = '  Estimated remaining time:'; call Message("(A$)")
 call PrintTime(TIME_unit_count*(float(numk)-float(ik)))
!
end subroutine Time_remaining
! ###################################################################
! 
!  subroutine PrintTime
!
!  Author: Marc De Graef
!  
!  Description: print the computation time
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
! ###################################################################
subroutine PrintTime(tm)

use local
use io

IMPLICIT NONE

integer(kind=irg)    :: days,hours,minutes,seconds
real(kind=sgl)       :: tm,secs

INTENT(IN)           :: tm

  secs = tm
  days = 0
  hours = 0
  minutes = 0
  if (secs.gt.86400.0) then
    days = int(secs)/86400
    secs = mod(secs,86400.0)
  end if
  if (secs.gt.3600.0) then
    hours = int(secs)/3600
    secs = mod(secs,3600.0)
  end if
  if (secs.gt.60.0) then
    minutes = int(secs)/60
    secs = mod(secs,60.0)
  end if
  seconds = int(secs)
  mess = ''; oi_int(1:4) = (/ days, hours, minutes, seconds /)
  call WriteInt(4,"(1x,I3,' d,',I3,' h,',I3,' m,',I3,' s'/)")
end subroutine PrintTime
! ###################################################################
! 
!  subroutine Time_stop
!
!  Author: Marc De Graef
!  
!  Description: end time recording
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/4/01  MDG 1.0 original
! ###################################################################
subroutine Time_stop(numk)

use local
use io

IMPLICIT NONE

integer(kind=irg)  :: numk

  call system_clock(TIME_newcount,TIME_count_rate,TIME_count_max)
  mess = '  Total computation time [s] '; call Message("(A$)")
  call PrintTime(float(TIME_newcount-TIME_count)/float(TIME_count_rate))
  mess = ' Time per step [s] ' 
  oi_real(1)=float(TIME_newcount-TIME_count)/float(TIME_count_rate)/float(numk)
  call WriteReal(1,"(F)")
end subroutine Time_stop

end module timing

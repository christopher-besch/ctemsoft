!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:  qg.f90                                                             !
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
! -Fortran-77 Source Code-
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "qg.f90"
!                                    created: 10/13/98 {9:29:46 AM} 
!                                 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description:  compute the complex extinction distance q_g
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/27/01 MDG 2.0 f90
! ###################################################################
program qg

use local
use io
use crystal
use files
use diffraction
use symmetry
use constants
use postscript
use dynamical

IMPLICIT NONE

integer(kind=irg)        :: ind(3),ans
real(kind=sgl)           :: preg

 progname = 'qg.f90'
 progdesc = 'Display potential coefficient values'
 call CTEMsoft
 
 call CrystalData
 call GetVoltage
 call CalcPositions('v')
 preg = 2.0 * sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18

 ans = 1
 do while (ans.eq.1)
  mess = 'Enter Miller indices :'; call Message("(/A)")
  call GetIndex(ind,'r')
  write (*,*) 'you selected ',ind
  call CalcUcg(ind)
  write (*,"(1x,'  h  k  l',1x,'   |Ug|    phase   |Ugp|   phase    xi_g   xi_gp   ratio   1/q_g')") 
  write (*,"(1x,3I3,1x,4F8.3,3F8.1,2F8.5)") rlp%hkl,rlp%Umod/preg,rlp%Vphase, &
                                                   rlp%Upmod/preg,rlp%Vpphase,rlp%xg,rlp%xgp,rlp%ar,rlp%qg
  write (*,"(1x,'Another one ? (1/0) : ',$)")
  read (*,"(I1)") ans
 end do 
end  program
       




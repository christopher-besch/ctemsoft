!----------------------------------------------------------------
! CTEMsoft                                                      !
! Copyright (c) 2001, Marc De Graef/Carnegie Mellon University  !
! All rights reserved.                                          !
!                                                               !
! The source code in this file is protected by a BSD license;   !
! for the complete text of the license see the file BSD.license.!
!                                                               !
!----------------------------------------------------------------
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "zap.f90"
!                                    created: 12/01/00 {9:29:46 AM} 
!                                last update: 12/01/00 {10:37:32 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: This program generates a multipage PostScript file
!               with zone axis patterns.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  12/11/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
program zap

use local
use crystalvars
use crystal
use symmetryvars
use symmetry
use graphics
use files
use postscript
use io
use diffraction

IMPLICIT NONE

integer(kind=irg)       :: inm

 progname = 'zap.f90'
 progdesc = 'Kinematical Zone Axis Diffraction Patterns'
 call CTEMsoft

 SG % SYM_reduce=.TRUE.
! read crystal information, microscope voltage, and camera length
 call CrystalData
 call GetVoltage
 mess = 'Camera length L  [mm, real] '; call GetReal(1)
 camlen = io_real(1)
! generate all atom positions in the fundamental unit cell
 call CalcPositions('v')
! open PostScript file
 call PS_openfile
! generate a set of zone axis patterns
 call DiffPage
! close Postscript file
 call PS_closefile
 call system(psviewer//PS%psname)
end program

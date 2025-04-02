! -*-Fortran-*-
! ###################################################################
! 
!  FILE: "crystal_variables.f90"
!                                    created: 1/5/99 {11:26:49 AM} 
!                                last update: 6/1/2001 {6:51:16 AM} 
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
!   4/ 5/00 MDG 1.2 modified TransCoor to include new mInvert
!   5/19/01 MDG 2.0 f90 version
! ###################################################################
!
!
! lattice parameters
! a             = a parameter [nm]
! b             = b parameter [nm]
! c             = c parameter [nm]
! alpha         = alpha angle [deg]
! beta          = beta  angle [deg]
! gamma         = gamma angle [deg]
! vol           = unit cell volume [nm^3]
!
! metric information
! dmt           = direct space metric tensor
! rmt           = reciprocal space metric tensor
! dsm           = direct space structure matrix
! rsm           = reciprocal space structure matrix
! krdel         = Kronecker delta (unit matrix)
!
!
! asymmetric unit contents
! ATOM_ntype    = actual number of occupied positions in asymmetric unit
! ATOM_type     = atomic number for each atom in asymmetric unit
! ATOM_pos      = fractional coordinates (x,y,z), occupation, Debye-Waller
!                 factor for each atom in asymmetric unit
! fname         = crystal structure file name
!

module crystalvars

use local

! the unitcell type contains everything that is needed for crystallographic
! computations, including the symmetry information for the space group
type unitcell
  real    :: a,b,c,alpha,beta,gamma,vol
  real    :: dmt(3,3),rmt(3,3),dsm(3,3),rsm(3,3),krdel(3,3)
  integer :: ATOM_type(maxpasym),ATOM_ntype,SYM_SGnum,xtal_system,SYM_SGset
  real    :: ATOM_pos(maxpasym,5)
  character(15) :: fname
  logical :: SYM_reduce,SYM_second,SYM_trigonal
end type

! this type is used to define an orientation relation, i.e., two parallel
! directions and two parallel planes
type orientation
  real    :: tA(3), tB(3), gA(3), gB(3)
end type

! cell is the generic unit cell variable used in all programs.  
type (unitcell) :: cell


end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:diffraction.f90                                                      !
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
!  FILE: "diffraction.f90"
!                                    created: 10/13/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: general diffraction routines
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
!   3/14/02 MDG 2.2 added CalcDynMat routine
! ###################################################################
module dynamical

use local

IMPLICIT NONE

! derived type definitions

type reflisttype  ! linked reflection list
  integer(kind=irg)          :: num, &  ! sequential number
                                hkl(3),&! Miller indices
				HOLZN,& ! belongs to this HOLZ layer
				famnum  ! family number
  character(1)               :: ForB    ! 'F' for full, 'B' for Bethe
  logical                    :: dbdiff  ! double diffraction reflection ?
  complex(kind=dbl)          :: Ucg,&   ! potential coefficient
                                amp     ! amplitude of beam
  real(kind=sgl)             :: sg      ! excitation error
  type(reflisttype),pointer  :: next    ! connection to next entry
end type reflisttype

type(reflisttype),pointer    :: reflist, & ! linked list of reflections
                                rltail, &  ! end of linked list
                                rltmpa,rltmpb ! temporary pointers

type result  ! linked wavevector list
  integer(kind=irg) :: i,j         ! image coordinates
  real(kind=sgl)    :: kt(3)       ! tangential component of wavevector
  real(kind=sgl)    :: kn          ! normal component
  real(kind=sgl)    :: k(3)        ! full wave vector
  type(result),pointer :: next     ! connection to next wave vector
end type result

type(result),pointer :: head, &    ! end of linked list
                        tmp, &     ! temporary pointer
                        tail       ! end of linked list

! allocatable arrays

complex(kind=dbl),allocatable :: W(:), &         ! eigenvalue vector for Bloch wave method
                                 CG(:,:), &      ! eigenvector matrix
                                 alpha(:), &     ! excitation amplitude vector
                                 DHWMz(:,:),&    ! Darwin-Howie-Whelan matrix
                                 DynMat(:,:), &  ! dynamical matrix
                                 phiz(:),Az(:,:) ! used for Taylor expansion of scattering matrix

! The parameters in gnode are computed by CalcUcg 
type gnode
  character(2)         :: method   ! computation method (WK = Weickenmeier-Kohl, DT = Doyle-Turner/Smith-Burge)
  logical              :: absorption ! is absorption included or not ?
  integer(kind=irg)    :: hkl(3)   ! Miller indices
  real(kind=sgl)       :: xg, &    ! extinction distance [nm]
                          xgp, &   ! absorption length [nm]
                          ar, &    ! aborption ratio
                          g, &     ! length of reciprocal lattice vectors [nm^-1]
                          Vmod,Vpmod, & ! modulus of Vg and Vgprime [V]
                          Umod,Upmod, & ! modulus of Ug and Ugprime [nm^-2]
                          Vphase,Vpphase ! phase factors of Vg and Vgprime [rad]
  complex(kind=sgl)    :: Ucg, &   ! scaled potential Fourier coefficient [nm^-2]
                          Vg, &    ! potential Fourier coefficient [V]
                          qg       ! interaction parameter for Darwin-Howie-Whelan equations [nm^-1]
end type gnode

type(gnode)            :: rlp      ! reciprocal lattice point

! other vectors needed for dynamical computations

real(kind=sgl)   :: DynWV(3), &       ! wave vector expressed in reciprocal frame
                    DynFN(3), &       ! Foil normal in reciprocal frame
		    DynUpz            ! U'_0 normal absorption parameter
integer(kind=irg):: DynNbeams      ! number of beams

end module dynamical


module doublediff

use local

! The following logical array is used to tag the potential double 
! diffraction spots in non-symmorphic space groups.

logical,allocatable  :: dbdiff(:)
logical              :: nonsymmorphic

end module doublediff
! 
! 
module multibeams

! these variables are passed back and forth between the main
! multi beam program and the various subroutines.  

use local

integer(kind=irg),parameter   :: numr = 500                    ! max number of families of reflections in zone
integer(kind=irg)             :: family(numr,48,3), numfam(numr) 
integer(kind=irg),allocatable :: idx(:)
real(kind=sgl)                :: glen(numr)                    ! length of g-vectors
real(kind=sgl),allocatable    :: gm(:),V(:,:)
logical,allocatable           :: al(:)                         ! array of allowed reflections


end module multibeams

! 
module diffraction

use local

!
! mAccvol       = microscope accelerating voltage  [V]
! mLambda       = electron wavelength [nm]
! mRelcor       = relativistic correction factor gamma [dimensionless]
! mSigma        = interaction constant [ ]
! mPsihat       = relativistic acceleration potential
! camlen        = diffraction camera length [mm]
! 

real(kind=sgl)           :: kzero(3),camlen
real(kind=dbl)           :: mAccvol,mLambda,mRelcor,mSigma,mPsihat
! atomic scattering factor parametrization (Doyle-Turner, Smith-Burge)
! used only if absorption is not taken into account;  otherwise
! the Weickenmeier-Kohl routine is used.
real(kind=sgl),parameter,private   :: scatfac(8,98) = reshape( (/ &
        0.202,0.244,0.082,0.000,0.30868,0.08544,0.01273,0.00000, &
        0.091,0.181,0.110,0.036,0.18183,0.06212,0.01803,0.00284, &
        1.611,1.246,0.326,0.099,1.07638,0.30480,0.04533,0.00495, &
        1.250,1.334,0.360,0.106,0.60804,0.18591,0.03653,0.00416, &
        0.945,1.312,0.419,0.116,0.46444,0.14178,0.03223,0.00377, &
        0.731,1.195,0.456,0.125,0.36995,0.11297,0.02814,0.00346, &
        0.572,1.043,0.465,0.131,0.28847,0.09054,0.02421,0.00317, &
        0.455,0.917,0.472,0.138,0.23780,0.07622,0.02144,0.00296, &
        0.387,0.811,0.475,0.146,0.20239,0.06609,0.01931,0.00279, &
        0.303,0.720,0.475,0.153,0.17640,0.05860,0.01762,0.00266, &
        2.241,1.333,0.907,0.286,1.08004,0.24505,0.03391,0.00435, &
        2.268,1.803,0.839,0.289,0.73670,0.20175,0.03013,0.00405, &
        2.276,2.428,0.858,0.317,0.72322,0.19773,0.03080,0.00408, &
        2.129,2.533,0.835,0.322,0.57775,0.16476,0.02880,0.00386, &
        1.888,2.469,0.805,0.320,0.44876,0.13538,0.02642,0.00361, &
        1.659,2.386,0.790,0.321,0.36650,0.11488,0.02469,0.00340, &
        1.452,2.292,0.787,0.322,0.30935,0.09980,0.02234,0.00323, &
        1.274,2.190,0.793,0.326,0.26682,0.08813,0.02219,0.00307, &
        3.951,2.545,1.980,0.482,1.37075,0.22402,0.04532,0.00434, &
        4.470,2.971,1.970,0.482,0.99523,0.22696,0.04195,0.00417, &
        3.966,2.917,1.925,0.480,0.88960,0.20606,0.03856,0.00399, &
        3.565,2.818,1.893,0.483,0.81982,0.19049,0.03590,0.00386, &
        3.245,2.698,1.860,0.486,0.76379,0.17726,0.03363,0.00374, &
        2.307,2.334,1.823,0.490,0.78405,0.15785,0.03157,0.00364, &
        2.747,2.456,1.792,0.498,0.67786,0.15674,0.03000,0.00357, &
        2.544,2.343,1.759,0.506,0.64424,0.14880,0.02854,0.00350, &
        2.367,2.236,1.724,0.515,0.61431,0.14180,0.02725,0.00344, &
        2.210,2.134,1.689,0.524,0.58727,0.13553,0.02609,0.00339, &
        1.579,1.820,1.658,0.532,0.62940,0.12453,0.02504,0.00333, &
        1.942,1.950,1.619,0.543,0.54162,0.12518,0.02416,0.00330, &
        2.321,2.486,1.688,0.599,0.65602,0.15458,0.02581,0.00351, &
        2.447,2.702,1.616,0.601,0.55893,0.14393,0.02446,0.00342, &
        2.399,2.790,1.529,0.594,0.45718,0.12817,0.02280,0.00328, &
        2.298,2.854,1.456,0.590,0.38830,0.11536,0.02146,0.00316, &
        2.166,2.904,1.395,0.589,0.33899,0.10497,0.02041,0.00307, &
        2.034,2.927,1.342,0.589,0.29999,0.09598,0.01952,0.00299, &
        4.776,3.859,2.234,0.868,1.40782,0.18991,0.03701,0.00419, &
        5.848,4.003,2.342,0.880,1.04972,0.19367,0.03737,0.00414, &
        4.129,3.012,1.179,0.000,0.27548,0.05088,0.00591,0.000,   &
        4.105,3.144,1.229,0.000,0.28492,0.05277,0.00601,0.000,   &
        4.237,3.105,1.234,0.000,0.27415,0.05074,0.00593,0.000,   &
        3.120,3.906,2.361,0.850,0.72464,0.14642,0.03237,0.00366, &
        4.318,3.270,1.287,0.000,0.28246,0.05148,0.00590,0.000,   &
        4.358,3.298,1.323,0.000,0.27881,0.05179,0.00594,0.000,   &
        4.431,3.343,1.345,0.000,0.27911,0.05153,0.00592,0.000,   &
        4.436,3.454,1.383,0.000,0.28670,0.05269,0.00595,0.000,   &
        2.036,3.272,2.511,0.837,0.61497,0.11824,0.02846,0.00327, &
        2.574,3.259,2.547,0.838,0.55675,0.11838,0.02784,0.00322, &
        3.153,3.557,2.818,0.884,0.66649,0.14449,0.02976,0.00335, &
        3.450,3.735,2.118,0.877,0.59104,0.14179,0.02855,0.00327, &
        3.564,3.844,2.687,0.864,0.50487,0.13316,0.02691,0.00316, &
        4.785,3.688,1.500,0.000,0.27999,0.05083,0.00581,0.000,   &
        3.473,4.060,2.522,0.840,0.39441,0.11816,0.02415,0.00298, &
        3.366,4.147,2.443,0.829,0.35509,0.11117,0.02294,0.00289, &
        6.062,5.986,3.303,1.096,1.55837,0.19695,0.03335,0.00379, &
        7.821,6.004,3.280,1.103,1.17657,0.18778,0.03263,0.00376, &
        4.940,3.968,1.663,0.000,0.28716,0.05245,0.00594,0.000,   &
        5.007,3.980,1.678,0.000,0.28283,0.05183,0.00589,0.000,   &
        5.085,4.043,1.684,0.000,0.28588,0.05143,0.00581,0.000,   &
        5.151,4.075,1.683,0.000,0.28304,0.05073,0.00571,0.000,   &
        5.201,4.094,1.719,0.000,0.28079,0.05081,0.00576,0.000,   &
        5.255,4.113,1.743,0.000,0.28016,0.05037,0.00577,0.000,   &
        6.267,4.844,3.202,1.200,1.00298,0.16066,0.02980,0.00367, &
        5.225,4.314,1.827,0.000,0.29158,0.05259,0.00586,0.000,   &
        5.272,4.347,1.844,0.000,0.29046,0.05226,0.00585,0.000,   &
        5.332,4.370,1.863,0.000,0.28888,0.05198,0.00581,0.000,   &
        5.376,4.403,1.884,0.000,0.28773,0.05174,0.00582,0.000,   &
        5.436,4.437,1.891,0.000,0.28655,0.05117,0.00577,0.000,   &
        5.441,4.510,1.956,0.000,0.29149,0.05264,0.00590,0.000,   &
        5.529,4.533,1.945,0.000,0.28927,0.05144,0.00578,0.000,   &
        5.553,4.580,1.969,0.000,0.28907,0.05160,0.00577,0.000,   &
        5.588,4.619,1.997,0.000,0.29001,0.05164,0.00579,0.000,   &
        5.659,4.630,2.014,0.000,0.28807,0.05114,0.00578,0.000,   &
        5.709,4.677,2.019,0.000,0.28782,0.05084,0.00572,0.000,   &
        5.695,4.740,2.064,0.000,0.28968,0.05156,0.00575,0.000,   &
        5.750,4.773,2.079,0.000,0.28933,0.05139,0.00573,0.000,   &
        5.754,4.851,2.096,0.000,0.29159,0.05152,0.00570,0.000,   &
        5.803,4.870,2.127,0.000,0.29016,0.05150,0.00572,0.000,   &
        2.388,4.226,2.689,1.255,0.42866,0.09743,0.02264,0.00307, &
        2.682,4.241,2.755,1.270,0.42822,0.09856,0.02295,0.00307, &
        5.932,4.972,2.195,0.000,0.29086,0.05126,0.00572,0.000,   &
        3.510,4.552,3.154,1.359,0.52914,0.11884,0.02571,0.00321, &
        3.841,4.679,3.192,1.363,0.50261,0.11999,0.02560,0.00318, &
        6.070,4.997,2.232,0.000,0.28075,0.04999,0.00563,0.000,   &
        6.133,5.031,2.239,0.000,0.28047,0.04957,0.00558,0.000,   &
        4.078,4.978,3.096,1.326,0.38406,0.11020,0.02355,0.00299, &
        6.201,5.121,2.275,0.000,0.28200,0.04954,0.00556,0.000,   &
        6.215,5.170,2.316,0.000,0.28382,0.05002,0.00562,0.000,   &
        6.278,5.195,2.321,0.000,0.28323,0.04949,0.00557,0.000,   &
        6.264,5.263,2.367,0.000,0.28651,0.05030,0.00563,0.000,   &
        6.306,5.303,2.386,0.000,0.28688,0.05026,0.00561,0.000,   &
        6.767,6.729,4.014,1.561,0.85951,0.15642,0.02936,0.00335, &
        6.323,5.414,2.453,0.000,0.29142,0.05096,0.00568,0.000,   &
        6.415,5.419,2.449,0.000,0.28836,0.05022,0.00561,0.000,   &
        6.378,5.495,2.495,0.000,0.29156,0.05102,0.00565,0.000,   &
        6.460,5.469,2.471,0.000,0.28396,0.04970,0.00554,0.000,   &
        6.502,5.478,2.510,0.000,0.28375,0.04975,0.00561,0.000,   &
        6.548,5.526,2.520,0.000,0.28461,0.04965,0.00557,0.000/), (/8,98/))

real(kind=dbl),allocatable    :: phir(:,:),phii(:,:),SMr(:,:,:),SMi(:,:,:)

real(kind=sgl),allocatable    :: Vg(:),rg(:),Vgsave(:)
integer(kind=irg),allocatable :: rfamily(:,:,:),rnumfam(:)

contains
! ###################################################################
!
!  subroutine GetVoltage
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 10/20/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: ask for acceleration voltage, and call CalcWaveLength
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine GetVoltage

use local
use io

IMPLICIT NONE

 mess ='Enter the microscope accelerating voltage [V] : '; call GetInt(1);
 call CalcWaveLength(dble(io_int(1)))
 
end subroutine
! ###################################################################
! 
!  subroutine CalcWaveLength
!
!  Author: Marc De Graef
!  
!  Description: computes the relativistic electron wavelength
!  These quantities are computed in double precision because of the 
!  wide range of magnitudes.  If a crystal structure has been defined
!  then the gamma*V_0 term is added to correct for refraction.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine CalcWaveLength(voltage,skip)

use local
use constants
use symmetry
use io
use dynamical

IMPLICIT NONE

real(kind=dbl)   :: voltage,temp1,temp2
real(kind=sgl)   :: Vmod,Vphase,Vpmod,Vpphase
integer(kind=irg):: hkl(3)
integer(kind=irg),OPTIONAL :: skip

intent(IN)       :: voltage,skip

! store voltage incommon block
 mAccvol = voltage
 temp1 = 1.0D+9*cPlanck/dsqrt(2.D0*cRestmass*cCharge)
 temp2 = cCharge*0.5D0*voltage/cRestmass/cLight**2
! relativistic correction factor (known as gamma)      
 mRelcor = 1.0D0+2.0D0*temp2
! relativistic acceleration voltage
 mPsihat = voltage*(1.D0+temp2)
! compute the electron wavelength in nm
! has a crystal structure been defined? If so, compute V_0 and
! add it to mPsihat (corrected by mRelcor)
 if (strucdef.eqv..TRUE.) then
  call CalcPositions('v')
! which scattering factors should be used ?
  if (present(skip)) then
   select case (skip) 
    case(1); rlp%method='DT'; 
    case(2); rlp%method='WK'; 
    case(3); rlp%method='WK'; rlp%absorption=.TRUE.
   end select
  else
   mess = ' The following scattering factor sets are available :'; call Message("(/A/)")
   mess = '  [1] Doyle-Turner/Smith-Burge (no absorption) '; call Message("(A)")
   mess = '  [2] Weickenmeier-Kohl (no absorption) '; call Message("(A)")
   mess = '  [3] Weickenmeier-Kohl (with absorption) '; call Message("(A/)")
   mess = 'Which set do you want to use [1/2/3] ? '; call GetInt(1)
   rlp%absorption = .FALSE.
   select case (io_int(1)) 
    case(1); rlp%method='DT'; 
    case(2); rlp%method='WK'; 
    case(3); rlp%method='WK'; rlp%absorption=.TRUE.
   end select
  end if
  hkl=(/0,0,0/)
  call CalcUcg(hkl) 
  mess = 'Mean inner potential [V] '; oi_real(1)=rlp%Vmod; call WriteReal(1,"(f8.4)")
  mPsihat = mPsihat + dble(rlp%Vmod)
  mess = ' Wavelength corrected for refraction'; call Message("(A)")
 endif
 mess = 'Relativistic correction factor [gamma]  '; oi_real(1)=sngl(mRelcor); call WriteReal(1,"(f14.6)")
 mess = 'Relativistic Accelerating Potential [V] '; oi_real(1)=sngl(mPsihat); call WriteReal(1,"(f14.6)")
 mLambda = temp1/dsqrt(mPsihat)
 mess = 'Electron Wavelength [nm]                '; oi_real(1)=sngl(mLambda); call WriteReal(1,"(f14.6)")
! interaction constant sigma
 mSigma = 2.D0*cPi*cRestmass*mRelcor*cCharge*mLambda
 mSigma = 1.0D-18*mSigma/cPlanck**2
 mess = 'Interaction constant [V nm]^(-1)        '; oi_real(1)=sngl(mSigma); call WriteReal(1,"(f14.6)")
end subroutine
!      
! ###################################################################
! 
!  function  CalcDiffAngle
!
!  Author: Marc De Graef
!  
!  Description: computes the Bragg scattering angle 2theta for a 
!               given set of Miller indices
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
real function CalcDiffAngle(h,k,l)

use local
use crystal

IMPLICIT NONE

integer(kind=irg)  :: h,k,l
real(kind=sgl)     :: ghkl, p(3)

! first compute the interplanar spacing
 p = float( (/h,k,l/) )
! then get the diffraction angle
 ghkl = CalcLength(p,'r')
 if (ghkl.ne.0.0) then
  ghkl = 2.0*asin(0.50*sngl(mLambda)*ghkl)
 end if
 CalcDiffAngle = ghkl
end function
! ###################################################################
! 
!  subroutine CalcUcg
!
!  Author: Marc De Graef
!  
!  Description: compute the complex Structure Factor for a given g
!               (CalcPositions MUST be called before this routine!)
!               includes extinction distance, absorption length, etc...
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine CalcUcg(hkl)

use local
use crystalvars
use crystal
use symmetry
use constants
use others
use dynamical

IMPLICIT NONE

integer(kind=irg)      :: hkl(3),j,absflg,m,jj,ii
real(kind=sgl)         :: s,twopi,arg,swk,dwwk,fp,pref,ul,pre,preg,sct,fs
complex(kind=sgl)      :: zz,ff,gg,sf,p1
complex(kind=sgl)      :: czero
logical                :: accflg, dwflg
character(2)           :: smb

intent(IN)             :: hkl

twopi=2.0*sngl(cPi)
czero = cmplx(0.0,0.0)
if (rlp%method.eq.'DT') then 
! compute the scattering parameter s^2=(g/2)^2
 if (hkl(1)**2+hkl(2)**2+hkl(3)**2.eq.0) then 
  s = 0.0
  rlp%g = 0.0
 else
  rlp%g = CalcLength(float(hkl),'r')
  s = (0.50*rlp%g)**2
 end if
! initialize the real and imaginary parts of the structure factor
 sf = czero
! compute the prefactor (this also scales from Angstrom to nm)
 pref = 0.04787801/cell % vol
! loop over all atoms in the asymmetric unit
 do m=1,cell % ATOM_ntype
! get the atomic scattering factor for this atom
  sct=0.0
  j=cell % ATOM_type(m)
  do ii=1,4
   sct=sct+scatfac(ii,j)*exp(-scatfac(ii+4,j)*s)
  end do
! scale and include Debye-Waller factor
! and site occupation parameter
  fs=pref*sct*exp(-cell % ATOM_pos(m,5)*s)*cell % ATOM_pos(m,4)
! loop over all atoms in the orbit
  do j=1,numat(m)
   arg=twopi*sum(hkl(1:3)*apos(m,j,1:3))
   sf = sf + fs*exp(cmplx(0.0,-arg))
  end do
 end do
! and fill in the entries of the rlp variable
 pre = 2.0*sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18
 rlp%hkl = hkl
 rlp%Vmod = cabs(sf)*mRelcor
 rlp%Vphase = atan2(aimag(sf),real(sf))
 rlp%Vpmod = 0.0
 rlp%Vpphase = 0.0
 if (rlp%Vmod.gt.0.0) then
  rlp%xg = 1.0/(pre*rlp%Vmod*mLambda)
 else
  rlp%xg = 1.0E+8
 end if
 rlp%xgp = 0.0
 rlp%ar = 0.0
 rlp%Vg = rlp%Vmod * exp(cmplx(0.0,rlp%Vphase))
 rlp%Ucg = pre*rlp%Vg
 rlp%qg = cmplx(1.0/rlp%xg,0.0)
else
! The Weickenmeier-Kohl (WK) subroutine works in Angstrom, and also 
! scales reciprocal space by a factor of 2*pi;  this scaling
! is accomplished by changing the g-value in nm^{-1} by a 
! scaling factor swk = 2*pi/10, to go from book units to WK units.
!
! A similar scaling must be performed on the Debye Waller factor;
! the book defines it as exp(-Bs^2), with s in [nm^2]; WK define
! it in A^2, and with a scaled reciprocal space.  The conversion
! factor dwwk = 100.0/8*pi^2
!
! To go from the standard B factor in [nm^2] to ul^2 in A^2,
! ul = sqrt(B*dwwk)
!
 swk = 0.1*twopi
 dwwk = 100.0/(8.0*cPi**2)
! compute the properly scaled length of g
 if (hkl(1)**2+hkl(2)**2+hkl(3)**2.eq.0) then 
  s = 0.0
  rlp%g = 0.0
 else
  rlp%g = CalcLength(float(hkl),'r')
  s = rlp%g*swk
 end if
! let fscatt perform the relativistic corrections for f_g and fprime_g
 accflg = .TRUE.
! include absorption ?
 absflg = 0
 if (rlp%absorption.eqv..TRUE.) then 
   absflg = 3  ! include phonon and core contributions
 end if
! always include Debye-Waller factor
 dwflg  = .TRUE.
! compute the scaling prefactors
! pref contains A to nm conversion, and divides by 4pi
 pref = 0.04787801/cell % vol/(4.0*cPi) 
! preg is used to go from V to U, remembering that gamma is already
! included in the output from fscatt
 preg = 2.0 * sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18
 pre = pref * preg
! initialize the real and imaginary parts of the structure factor
 ff=czero
 gg=czero
! loop over all atoms in the asymmetric unit
 do m=1,cell % ATOM_ntype
! get the atomic scattering factor for this atom
! scale and include Debye-Waller factor and site occupation parameter
  ul = sqrt(cell % ATOM_pos(m,5)*dwwk)
  j = cell % ATOM_type(m)
  sf = fscatt(s,ul,j,smb,sngl(mAccvol)/1000.0,absflg,accflg,dwflg)*cell%ATOM_pos(m,4)
! loop over all atoms in the orbit
  p1 = czero
  do j=1,numat(m)
   arg=twopi*sum(float(hkl(1:3))*apos(m,j,1:3))
   p1 = p1 + exp(cmplx(0.0,-arg))
  end do
  ff = ff + p1*real(sf)
  gg = gg + p1*aimag(sf)
 end do
!
! fill in the entries of the rlp variable
 rlp%hkl = hkl
! these are the modulus and phase of the real part of Vg
 rlp%Vmod = pref * cabs(ff)
 rlp%Vphase = atan2(aimag(ff),real(ff))
! modulus of U_g
 rlp%Umod = preg*rlp%Vmod
! if absorption is included, also compute the imaginary part of Vg, i.e., Vprime_g
 if (rlp%absorption.eqv..TRUE.) then 
  rlp%Vpmod = pref * cabs(gg)
  rlp%Vpphase = atan2(aimag(gg),real(gg))
! modulus of Uprime_g
  rlp%Upmod = preg*rlp%Vpmod
! complex Ucg = U_g + i Uprime_g = U_g,r-Uprime_g,i + i(U_g,i+Uprime_g,r)
  rlp%Ucg = pre * cmplx(real(ff)-aimag(gg),aimag(ff)+real(gg))
 else ! set absorption parameters to zero
  rlp%Vpmod = 0.0
  rlp%Vpphase = 0.0
! Ucg = U_g (complex number)
  rlp%Ucg = pre * ff
 end if
! complex Vg 
 rlp%Vg = rlp%Ucg/preg
 if (abs(rlp%Umod).gt.0.0) then 
  rlp%xg = 1.0/abs(rlp%Umod)/mLambda
 else
  rlp%xg = 1.0E+8
 end if 
 if (abs(rlp%Upmod).gt.0.0) then 
  rlp%xgp = 1.0/abs(rlp%Upmod)/mLambda
 else
  rlp%xgp = 1.0E+8
 end if 
 if (rlp%absorption.eqv..TRUE.) then 
  rlp%ar = rlp%xgp/rlp%xg
  arg = rlp%Vpphase-rlp%Vphase
  rlp%qg = cmplx(1.0/rlp%xg-sin(arg)/rlp%xgp,cos(arg)/rlp%xgp)
 else
  rlp%ar = 0.0
  rlp%qg = cmplx(1.0/rlp%xg,0.0)
 end if
end if
end subroutine
! ###################################################################
! 
!  subroutine Printrlp
!
!  Author: Marc De Graef
!  
!  Description: output the contents of the rlp structure
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  12/08/01 MDG 1.0 original
! ###################################################################
subroutine Printrlp(first)

use local
use constants
use dynamical

IMPLICIT NONE

logical,optional,intent(IN)  :: first

if (present(first)) then
 write (*,"(//'scattering factors : ',$)")
 if (rlp%method.eq.'WK') then 
  if (rlp%absorption.eqv..TRUE.) then 
   write (*,"('Weickenmeier-Kohl (with absorption)'/)")
  else
   write (*,"('Weickenmeier-Kohl'/)")
  end if
 else
   write (*,"('Doyle-Turner/Smith-Burge'/)")
 end if
 if (rlp%absorption.eqv..TRUE.) then
   write (*,"(1x,'  h  k  l',1x,'  |g|    Ucg_r  Ucg_i',1x,'   |Ug|    phase   |Ugp|   phase    xi_g   xi_gp   ratio   1/q_g')") 
 else
   write (*,"(1x,'  h  k  l',1x,'  |g|    Ucg_r',1x,'   |Ug|    phase    xi_g   1/q_g')") 
 end if
end if
if (rlp%absorption.eqv..TRUE.) then
 write (*,"(1x,3I3,1x,F7.4,2F7.3,1x,4F8.3,3F8.1,2F8.3)") rlp%hkl,rlp%g,rlp%Ucg,rlp%Umod,rlp%Vphase*180.0/cPi, &
                                                   rlp%Upmod,rlp%Vpphase*180.0/cPi,rlp%xg,rlp%xgp,rlp%ar,rlp%qg
else
 write (*,"(1x,3I3,1x,F7.4,F7.3,1x,2F8.3,F8.1,2F8.3)") rlp%hkl,rlp%g,real(rlp%Ucg),rlp%Umod,rlp%Vphase*180.0/cPi,rlp%xg,rlp%qg
end if

end subroutine
! ###################################################################
! 
!  function Calcsg
!
!  Author: Marc De Graef
!  
!  Description: compute excitation error sg for a given reflection
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
real function Calcsg(gg,kk,FN)

use local
use crystal

IMPLICIT NONE

real(kind=sgl) ::  gg(3),kk(3),kpg(3),tkpg(3),xnom,xden,q1,q2,FN(3)

intent(IN)     :: gg,kk,FN

 kpg=kk+gg
 tkpg=2.0*kk+gg
! use equation of Ewald sphere
 xnom = -CalcDot(gg,tkpg,'r')
! 2|k0+g|cos(alpha) = 2(k0+g).Foilnormal
 q1 = CalcLength(kpg,'r')
 q2 = CalcAngle(kpg,FN,'r')
 xden = 2.0*q1*cos(q2)
 Calcsg = xnom/xden
end function
! ###################################################################
! 
!  subroutine CalcDynMat
!
!  Author: Marc De Graef
!  
!  Description: compute the dynamical matrix for either Bloch wave
!               or Darwin-Howie-Whelan approach. Bethe potentials
!               are not yet implemented in this routine, but should
!               be at some point in the future...
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  3/14/02 MDG 1.0 original
! ###################################################################
subroutine CalcDynMat(calcmode)

use local
use constants
use error
use crystal
use dynamical

IMPLICIT NONE

character(5)       :: calcmode  ! D-H-W, BLOCH, DIAGH, DIAGB, or BETHE
integer(kind=irg)  :: istat,ir,ic
real(kind=sgl)     :: glen,exer
complex(kind=dbl)  :: czero,pre

intent(IN)         :: calcmode

! has reflist been allocated ?
if (.not.associated(reflist)) then 
  call FatalError('CalcDynMat: reflection list has not been allocated',' ')
end if
! initialize some parameters
czero = cmplx(0.0,0.0,dbl)
pre = cmplx(0.0,cPi,dbl)
! allocate DynMat if it hasn't already been allocated
if (.not.allocated(DynMat)) then
  allocate(DynMat(DynNbeams,DynNbeams),stat=istat)
  DynMat = czero
! get the absorption coefficient
  call CalcUcg((/0,0,0/))
  DynUpz = rlp%Vpmod
end if
! are we supposed to fill the off-diagonal part ?
 if ((calcmode.eq.'D-H-W').or.(calcmode.eq.'BLOCH')) then
  rltmpa => reflist%next    ! point to the front of the list
! ir is the row index
  do ir=1,DynNbeams
   rltmpb => reflist%next   ! point to the front of the list
! ic is the column index
   do ic=1,DynNbeams
    if (ic.ne.ir) then  ! exclude the diagonal
! compute Fourier coefficient of electrostatic lattice potential 
     call CalcUcg(rltmpa%hkl - rltmpb%hkl)
     if (calcmode.eq.'D-H-W') then
      DynMat(ir,ic) = pre*rlp%qg
     else
      DynMat(ir,ic) = rlp%Ucg
     end if
    end if
    rltmpb => rltmpb%next  ! move to next column-entry
   end do
   rltmpa => rltmpa%next   ! move to next row-entry
  end do
 end if
! or the diagonal part ?
 if ((calcmode.eq.'DIAGH').or.(calcmode.eq.'DIAGB')) then
  rltmpa => reflist%next   ! point to the front of the list
! ir is the row index
  do ir=1,DynNbeams
   glen = CalcLength(float(rltmpa%hkl),'r')
   if (glen.eq.0.0) then
    DynMat(ir,ir) = czero
   else  ! compute the excitation error
    exer = Calcsg(float(rltmpa%hkl),DynWV,DynFN)
    rltmpa%sg = exer
    if (calcmode.eq.'DIAGH') then  !
     DynMat(ir,ir) = cmplx(0.0,2.D0*cPi*exer,dbl)
    else
     DynMat(ir,ir) = cmplx(2.D0*exer/mLambda,DynUpz,dbl)
    end if
   endif
   rltmpa => rltmpa%next   ! move to next row-entry
  end do
 end if
! or are we making use of Bethe potentials ?
! This portion is yet to be implemented
end subroutine CalcDynMat
! ###################################################################
! 
!  subroutine MakeRefList
!
!  Author: Marc De Graef
!  
!  Description: allocate and initialize the linked reflection list
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  3/14/02 MDG 1.0 original
! ###################################################################
subroutine MakeRefList

use local
use error
use dynamical

IMPLICIT NONE

integer(kind=irg)  :: istat

! create it if it does not already exist
if (.not.associated(reflist)) then
  DynNbeams = 0
  allocate(reflist,stat=istat)
  if (istat.ne.0) call FatalError('MakeRefList: unable to allocate pointer',' ')
  rltail => reflist                 ! tail points to new value
  nullify(rltail%next)              ! nullify next in new value
end if
end subroutine MakeRefList
! ###################################################################
! 
!  subroutine AddReflection
!
!  Author: Marc De Graef
!  
!  Description: add a reflection to the linked reflection list
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  3/14/02 MDG 1.0 original
! ###################################################################
subroutine AddReflection(hkl)

use local
use error
use dynamical

IMPLICIT NONE

integer(kind=irg)  :: istat,hkl(3)

intent(IN)         :: hkl

! create reflist if it does not already exist
 if (.not.associated(reflist)) call MakeRefList
! create a new entry
 allocate(rltail%next,stat=istat)  ! allocate new value
 if (istat.ne.0) call FatalError('AddReflection: unable to add new reflection',' ')
 rltail => rltail%next             ! tail points to new value
 nullify(rltail%next)              ! nullify next in new value
 DynNbeams = DynNbeams + 1         ! update reflection counter
 rltail%num = DynNbeams            ! store reflection number
 rltail%hkl = hkl                  ! store Miller indicees
 call CalcUcg(hkl)                 ! compute potential Fourier coefficient
 rltail%Ucg = rlp%Ucg              ! and store it in the list
end subroutine AddReflection
! ###################################################################
! 
!  subroutine TBCalcSM
!
!  Author: Marc De Graef
!  
!  Description: compute crystal scattering matrix (real and imaginary
!  part) for a given set of parameters. Optimized to minimize the 
!  total number of function evaluations and multiplications/divisions
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine TBCalcSM(Ar,Ai,sg,z,xig,xigp,xizero,betag)

use local
use constants

IMPLICIT NONE
        
real(kind=sgl)  :: Ar(2,2), Ai(2,2), sg, z, pr, pi, cs, ss, ch, sh, q, q1, q2, sgs, &
                   sr, si, o , p, sb, cb, e, xig, xigp, xizero, betag, r, sq, xigi,xigpi
     
intent(IN)      :: sg,z,xig,xigp,xizero,betag
intent(OUT)     :: Ar,Ai

! setup auxiliary variables 
 xigi = 1.00/xig
 xigpi = 1.00/xigp
! sigma squared
 cb = cos(betag)
 sb = sin(betag)
 q = sg**2+xigi**2-xigpi**2
 r = cb*xigi*xigpi
 sq = sqrt(q**2+4.0*r**2)
! real part of sigma
 sr = sqrt(0.5*(q+sq))
! imaginary part of sigma
 si = r/sr
 sq = 1.0/sq
! s_g divided by sigma squared
 sgs = sg*sq
! arguments of trigonometric and hyperbolic functions
 e = cPi*z
 pr = e*sr
 pi = e*si
 e = exp(-e/xizero)
! trigonometric and hyperbolic functions
 cs = cos(pr)
 ss = sin(pr)
 ch = cosh(pi)
 sh = sinh(pi)
 o = ss*ch
 p = cs*sh
! transmitted amplitude T including normal absorption
 q = e*sgs
 q1 = q*(si*o-sr*p)
 q2 = q*(sr*o+si*p)
 Ar(1,1) = e*cs*ch
 Ai(1,1) = -e*ss*sh
 Ar(2,2) = Ar(1,1)+q1
 Ai(2,2) = Ai(1,1)+q2
 Ar(1,1) = Ar(1,1)-q1
 Ai(1,1) = Ai(1,1)-q2
! scattered amplitude S including normal absorption
 q1 = e*sq*(si*xigi-(sr*cb+si*sb)*xigpi)
 q2 = e*sq*(sr*xigi-(sr*sb-si*cb)*xigpi)
 Ar(1,2) = q1*o-q2*p
 Ai(1,2) = q2*o+q1*p
 Ar(2,1) = Ar(1,2)
 Ai(2,1) = Ai(1,2)
end subroutine
! ###################################################################
! 
!  subroutine TBCalcInten
!
!  Author: Marc De Graef
!  
!  Description: compute transmitted and scattered intensities for
!  the perfect crystal case.   This routine does not make use of 
!  complex number arithmetic, but instead uses the analytical 
!  expressions for the two-beam intensities derived in section 
!  6.3.3.4 on page 356-357.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine TBCalcInten(It,Is,sg,z,xig,xigp,xizero,betag)

use local
use constants

IMPLICIT NONE

real(kind=sgl) :: It, Is, sg, z, xig, xigp, xizero, betag, q, r, sq, qgsi, e, &
                  sr, si, cp, ch, pr, pi, xigi, xigpi, sgs
     
intent(IN)     :: sg,z,xig,xigp,xizero,betag
intent(OUT)    :: It,Is

! setup auxiliary variables 
 xigi = 1.0/xig
 xigpi = 1.0/xigp
! sigma squared
 q = sg**2+xigi**2-xigpi**2
 r = cos(betag)*xigi*xigpi
 sq = sqrt(q**2+4.0*r**2)
! real part of sigma
 sr = sqrt(0.5*(q+sq))
! imaginary part of sigma
 si = r/sr
 sq = 1.0/sq
! reciprocal of q_g squared
 qgsi = xigi**2+xigpi**2-2.0*sin(betag)*xigi*xigpi
! s_g squared divided by sigma squared
 sgs = sg**2*sq
! arguments of trigonometric and hyperbolic functions
 e = 2.0*cPi*z
 pr = e*sr
 pi = e*si
 e = exp(-e/xizero)
! trigonometric functions 
 cp = cos(pr)
 ch = cosh(pi)
! transmitted intensity It
 It = 0.5*ch*(1.0+sgs)+sg*sq*(sr*sinh(pi)-si*sin(pr))+0.5*cp*(1.0-sgs)
 It = It*e
! scattered intensity Is
 Is = 0.5*qgsi*sq*e*(ch-cp)
end subroutine
! ###################################################################
! 
!  subroutine TBCalcdz 
!
!  Author: Marc De Graef
!  
!  Description: multiply the column vector with the Scattering 
!  Matrix.                  
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine TBCalcdz(im,nbm)
 
use local

IMPLICIT NONE

real(kind=dbl),allocatable :: p(:),q(:)
integer(kind=irg)          :: im,nbm,k,i,j

intent(IN)                 :: im,nbm

 allocate(p(nbm))
 allocate(q(nbm))
 do k=1,im
  do i=1,nbm
   p(i)=0.D0
   q(i)=0.D0
   do j=1,nbm
    p(i)=p(i)+SMr(k,i,j)*phir(k,j)-SMi(k,i,j)*phii(k,j)
    q(i)=q(i)+SMi(k,i,j)*phir(k,j)+SMr(k,i,j)*phii(k,j)
   end do
  end do
  do i=1,nbm
   phir(k,i)=p(i)
   phii(k,i)=q(i)
  end do
 end do
 deallocate(p)
 deallocate(q)
end subroutine
! ###################################################################
!
!  subroutine DiffPage
!
!  Author: Marc De Graef
!
!  Description: draw zone axis electron diffraction patterns
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine DiffPage

use local
use postscript
use crystal
use crystalvars
use symmetry
use symmetryvars
use math
use io
use constants
use doublediff
use dynamical

IMPLICIT NONE


integer(kind=irg),parameter     :: inm = 5
character(1)                    :: list(256)
logical                         :: first,np,ppat,indi,a
logical,allocatable             :: z(:,:,:),zr(:,:,:)
integer(kind=irg)               :: i,j,h,k,l,m,totfam,hh,ll,fmax,inmhkl(3),ricnt,icnt,ind(3),uu,vv,ww,slect(256), &
                                   js,ii,num,hc,hhcc,iinm,dpcnt,imo,ih,ik,il,ier,iref
integer(kind=irg),allocatable   :: idx(:)
integer(kind=irg),allocatable   :: family(:,:),numfam(:)
real(kind=sgl)                  :: twopi,ggl,g(3),Vmod,Vphase,Vpmod,Vpphase,igl,Vmax,laL,gmax,RR,thr
real(kind=sgl),allocatable      :: gg(:)
real(kind=sgl),parameter        :: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/), &
                                   eps = 1.0E-3

 SG % SYM_reduce=.TRUE.
 thr = 1.E-4 
 twopi = 2.0*cPi
 Vmax = 0.0
! determine the families of reciprocal lattice points
! in a region of reciprocal space.
! gmax is the radius of the sphere whose intersection with the 
! back focal plane is the circle on the output zone axis patterns
 laL = sngl(mLambda) * camlen
 mess = 'wavelength [nm] = '; oi_real(1) = sngl(mLambda); call WriteReal(1,"(F10.6)")
 mess = ' L         [mm] = '; oi_real(1) = camlen; call WriteReal(1,"(f10.2)")
 mess = 'camera length lambda*L [mm nm] = '; oi_real(1) = laL; call WriteReal(1,"(f10.5)")
 RR = 1.375 * 25.4
 gmax = RR/laL
! set the index boundaries
 do i=1,3
  inmhkl(i) = 2*int(gmax/sqrt(cell % rmt(i,i)))
 end do
 hc = maxval(inmhkl)
! allocate all the necessary arrays
 allocate(zr(-hc:hc,-hc:hc,-hc:hc))
 allocate(z(-inm:inm,-inm:inm,-inm:inm))
 hhcc = (2*hc+1)**3
 allocate(Vg(hhcc))
 allocate(Vgsave(hhcc))
 allocate(rfamily(hhcc,48,3))
 allocate(rnumfam(hhcc))
 allocate(rg(hhcc))
 iinm = (2*inm+1)**3
 allocate(family(iinm,3))
 allocate(numfam(iinm))
! if this is a non-symmorphic space group, then also
! allocate the dbdiff array to tag potential double 
! diffraction reflections
 nonsymmorphic = (minval(abs(SGsym - cell % SYM_SGnum)).ne.0) 
 if (nonsymmorphic) then
   allocate(dbdiff(hhcc))
   dbdiff(1:hhcc) = .FALSE.
 endif
! and initialize the ones that need to be initialized
 zr(-hc:hc,-hc:hc,-hc:hc) = .FALSE.
 z(-inm:inm,-inm:inm,-inm:inm) = .FALSE.
! here we go:
 first = .TRUE.
 icnt = 1
 totfam=0
 do h=-inmhkl(1),inmhkl(1)
  ind(1)=h
  do k=-inmhkl(2),inmhkl(2)
   ind(2)=k
   do l=-inmhkl(3),inmhkl(3)
    ind(3)=l
! make sure we have not already done this one in another family
    if (.not.zr(h,k,l)) then
! check the length, to make sure it lies within the sphere gmax
     ggl=CalcLength(float(ind),'r')
! if it is larger than gmax, then compute the entire family
     if (ggl.ge.gmax) then
      call CalcFamily(ind,num,'r')
! and label the family members in the zr array so that we 
! do not include those points later on in the loop
      do i=1,num
       zr(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
      end do
     else 
! if the length is smaller than gmax, then compute the 
! Fourier coefficient Vg and determine the entire family
! [recall that all members in a family have the same Vg]
! Do this only for those reflections that are allowed by
! the lattice centering !
      a = IsGAllowed((/h,k,l/))
      if (a) then
       call CalcUcg(ind)
! check for nonsymmorphic systematic absences
       if ((nonsymmorphic).and.(rlp%Vmod.lt.eps)) then
        mess = 'potential double diffraction family :'
        oi_int(1) = h
        oi_int(2) = k
        oi_int(3) = l
        call WriteInt(3,"('{',I3,I3,I3,'}')")
        dbdiff(icnt) = .TRUE.
        rlp%Vmod = 0.0
       endif
! compute the entire family
       call CalcFamily(ind,num,'r')
       rg(icnt)=ggl
! copy family in array
       do i=1,num
        rfamily(icnt,i,1:3)=itmp(i,1:3)
        zr(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
       end do
! and take the modulus squared for the intensity
       Vg(icnt)=rlp%Vmod**2
       Vgsave(icnt)=Vg(icnt)
! update the maximum value 
       Vmax = max(Vg(icnt),Vmax)
      else
! remove the equivalent systematic absences
       call CalcFamily(ind,num,'r')
       rg(icnt)=ggl
       do i=1,num
        rfamily(icnt,i,1:3)=itmp(i,1:3)
        zr(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
       end do
! and put the intensity to a negative value
! that way we will know whether or not to draw them
       Vg(icnt)=-100.0
       Vgsave(icnt)=Vg(icnt)
      end if
! and increment the family counter
      rnumfam(icnt)=num
      totfam=totfam+num-1
      icnt=icnt+1
     end if 
    end if 
   end do
  end do
 end do
 icnt=icnt-1
! normalize potential coefficients to largest one
! and scale in a non-linear way to mimic density on 
! an electron micrograph [Gonzalez & Windtz]
! Use the where operator to avoid the negative intensities
 mess = 'logarithmic[0] or exponential[1] intensity scale'; call GetInt(1)
 ll = io_int(1)
 if (ll.eq.0) then
  where(Vg.gt.0.0) Vg=0.05*alog(1.0+0.1*Vg)
 else
  where(Vg.gt.0.0) Vg=0.05*(Vg/Vmax)**0.2
 end if
 ricnt=icnt
! determine families of directions
 first = .TRUE.
 icnt = 1
 totfam=0
 do uu=-inm,inm
  do vv=-inm,inm
   do ww=-inm,inm
    if ((uu**2+vv**2+ww**2).ne.0) then
! make sure we have not already done this one in another family
     ind= (/ -uu, -vv, -ww /)
     call IndexReduce(ind)
     if (.not.z(ind(1),ind(2),ind(3))) then
! determine the family <uvw>
      call CalcFamily(ind,num,'d')
! and keep only one family member, namely the one with the
! largest sum of the three integers, i.e. u+v+w
! [this is a simple way to get mostly positive indices as
! the zone axis indices]
      js = -100
      ii = 0
      do i=1,num
       hh = itmp(i,1)+itmp(i,2)+itmp(i,3)
       if (hh.gt.js) then 
        ii = i
        js = hh
       end if
! then remove the multiples of those direction indices from list 
       do m=-inm,inm
        ih=itmp(i,1)*m
        ik=itmp(i,2)*m
        il=itmp(i,3)*m
        if (((abs(ih).le.inm).and.(abs(ik).le.inm)).and.(abs(il).le.inm)) then 
         z(ih,ik,il)=.TRUE.
        end if
       end do
      end do
      family(icnt,1:3)=itmp(ii,1:3)
! increment family counter
      numfam(icnt)=num
      totfam=totfam+num-1
      icnt=icnt+1
     end if
    end if
   end do
  end do
 end do
 icnt=icnt-1
 mess = '->Total number of direction families = '; oi_int(1)=icnt; call WriteInt(1,"(I6)")
! compute length of direction vectors and rank by increasing length
 allocate(idx(icnt))
 allocate(gg(icnt))
 gg(1:icnt) = 0.0
 do k=1,icnt
  g(1:3)=float(family(k,1:3))
  gg(k)=CalcLength(g,'d')
 end do
! rank by increasing value of gg (use SLATEC routine)
 call SPSORT(gg,icnt,idx,1,ier)
! ask for number to be included in output
 mess = 'List of available zone axis patterns'; call Message("(A)")
 do i=1,icnt
  j=idx(i)
  oi_int(1)=i
  do k=1,3
   oi_int(k+1) = family(j,k)
  end do
  if (mod(i,4).eq.0) then 
   call WriteInt(4,"(I3,' [',3I3,'];')")
  else
   call WriteInt(4,"(I3,' [',3I3,'];',$)")
  endif
 end do
 mess = 'Enter selection (e.g. 4,10-20,) '; call Message("(//,A)")
 mess = '[Include 0 to also draw a powder pattern] '; call Message("(A)")
 read (*,fmt="(256A)") list
 call studylist(list,slect,fmax,ppat)
 if (nonsymmorphic) then
  mess = 'Potential double diffraction reflections will be indicated by open squares.'; call Message("(A,/)")
 end if
 mess = 'No indices (0), labels (1), extinctions (2), labels + extinctions (3): '; call GetInt(1)
 iref = io_int(1)
! and create output in 2 columns, 3 rows 
 do i=1,fmax  
  dpcnt=dpcnt+1
  j=idx(slect(i))
  imo = mod(i-1,6)
  if (imo.eq.0) then 
   np=.TRUE.
  else
   np=.FALSE.
  endif
  if (i.eq.1) then
   first=.TRUE.
  else
   first=.FALSE.
  endif
  if (slect(i).eq.0) then
   mess = 'Creating Powder Pattern '; call Message("(A)")
   call DumpPP(xoff(imo),yoff(imo),np,laL,ricnt)
   ppat=.FALSE.
  else
   mess = 'Creating ZAP ';
   do k=1,3
     oi_int(k) = family(j,k)
   end do
   call WriteInt(3,"('[',3i3,'] : ',$)")
   call DumpZAP(xoff(imo),yoff(imo),family(j,1),family(j,2),family(j,3),numfam(j),np,first,iref,laL,ricnt)
  endif
 end do
! and clean up all variables
 deallocate(zr)
 deallocate(z)
 deallocate(Vg)
 deallocate(Vgsave)
 deallocate(rfamily)
 deallocate(rnumfam)
 deallocate(rg)
 deallocate(family)
 deallocate(numfam)
 deallocate(idx)
 deallocate(gg)
 if (nonsymmorphic) deallocate(dbdiff)
end subroutine
! ###################################################################
!
!  subroutine DumpZAP
!
!  Author: Marc De Graef
!
!  Description: draw a single zone axis diffraction pattern
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine DumpZAP(xo,yo,u,v,w,p,np,first,indi,laL,icnt)

use local
use io
use postscript
use crystal
use crystalvars
use doublediff
use error

IMPLICIT NONE

! nref is the anticipated maximum number of reflections per pattern
integer(kind=irg),parameter  :: nref = 2000
integer(kind=irg)            :: u,v,w,p,dp,indi,i,j,jcnt,ui,vi,wi,pp,icnt,locg(nref,3),ier
integer(kind=irg),allocatable:: idx(:)
real(kind=sgl)               :: xo,yo,sc,laL,gmax,leng(nref),PX,PY,qx,qy,locv(nref),locvsave(nref),t(3),c(3),gg(3),gx(3),gy(3)
real(kind=sgl),allocatable   :: lng(:)
real(kind=sgl),parameter     :: le=3.25,he=2.9375,thr=1.0E-4
logical                      :: np,first,dbd(nref)

! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
 if (np) then
  call PS_newpage(.FALSE.,'Kinematical Zone Axis Patterns')
  call PS_text(5.25,-0.05,'scale bar in reciprocal nm')
  gmax = laL
  call PS_textvar(5.25,PS % psfigheight+0.02,'Camera Constant [nm mm]',gmax)
  call PS_setfont(PSfonts(2),0.15)
  call PS_text(-0.25,PS % psfigheight+0.02,'Structure File : '//cell % fname)
 end if
! draw frame and related stuff
 call PS_setlinewidth(0.012)
 call PS_balloon(xo,yo,le,he,0.0312)
 call PS_setlinewidth(0.001)
 PX = xo+1.8125
 PY = yo+1.0+15.0/32.0
 call PS_circle(PX,PY,1.375)
! zone axis
 call PS_setfont(PSfonts(2),0.12)
 call PS_text(xo+0.05,yo+he-0.15,'Zone axis ')
 ui=u
 vi=v
 wi=w
 call PrintIndices('d',ui,vi,wi,xo+0.6,yo+he-0.15)
! multiplicity
 call PS_setfont(PSfonts(2),0.12)
 pp=p
 call PS_textint(xo+0.05,yo+he-0.30,'Multiplicity ',pp)
! scale bar (sc is the conversion factor from nm-1 to inches)
 sc = laL/25.4
 call PS_setlinewidth(0.020)
 call PS_line(xo+0.05,yo+0.06,xo+0.05+5.0*sc,yo+0.06)
 call PS_setfont(PSfonts(2),0.15)
 call PS_text(xo+0.05+2.5*sc,yo+0.10,'5 ')
! select all reflections belonging to this zone axis pattern
 leng(1:nref)=0.0
 jcnt=0
 do i=1,icnt
  do j=1,rnumfam(i)
   dp=u*rfamily(i,j,1)+v*rfamily(i,j,2)+w*rfamily(i,j,3)
   if (dp.eq.0) then
    jcnt=jcnt+1
    if (jcnt.gt.nref) call Errormess(7)
    locg(jcnt,1:3)=rfamily(i,j,1:3)
    leng(jcnt)=rg(i)
    locv(jcnt)=Vg(i)
    locvsave(jcnt)=Vgsave(i)
    dbd(jcnt)=.FALSE.
! take care of potential double diffraction reflections
    if ((nonsymmorphic).and.(dbdiff(i))) dbd(jcnt) = .TRUE.
   end if
  end do
 end do
! rank them by length (use SLATEC routine)
 allocate(idx(jcnt))
 allocate(lng(jcnt))
 lng(1:jcnt) = leng(1:jcnt)
 call SPSORT(lng,jcnt,idx,1,ier)
 mess = ' Number of reflections : '; oi_int(1)=jcnt; call WriteInt(1,"(I5)")
! normalize the zone axis in cartesian components; this is the z-axis
 t(1)=-float(u)
 t(2)=-float(v)
 t(3)=-float(w)
 call TransSpace(t,c,'d','c')
 call NormVec(c,'c')
! take the first reflection in the list and make that the x-axis
! skip the zero reflection !!
 j=idx(2)
! normalize the first reciprocal vector in cartesian components
! this will be the x-axis of the diffraction pattern
 gg(1:3)=float(locg(j,1:3))
 call TransSpace(gg,gx,'r','c')
 call NormVec(gx,'c')
! then get the cross product between t and g; this is the y-axis
 call CalcCross(c,gx,gy,'c','c',0)
! plot origin of reciprocal space 
 call PS_filledcircle(PX,PY,0.05,0.0)
! then plot the remaining reflections
 do i=1,jcnt
  j=idx(i)
  gg(1:3)=float(locg(j,1:3))
  call TransSpace(gg,c,'r','c')
  qx=PX-CalcDot(c,gx,'c')*sc
  qy=PY+CalcDot(c,gy,'c')*sc
! first check for systematic absence due to lattice centering
  if ((locvsave(j).eq.-100.0).and.(indi.ge.2)) then
    call PS_cross(qx,qy,0.03,0.001)
  end if
! could it be a double diffraction spot ?
  if ((nonsymmorphic).and.(dbd(j))) call PS_square(qx,qy,0.04)
! is it a regular reflection ?
  if (locv(j).ge.thr) then
   call PS_filledcircle(qx,qy,locv(j),0.0)
   if ((indi.eq.1).or.(indi.eq.3)) then 
    call Printhkl(qx,qy,locg(j,1),locg(j,2),locg(j,3))
   end if
  end if
 end do
 deallocate(idx)
 deallocate(lng)
end subroutine
! ###################################################################
!
!  subroutine DumpPP
!
!  Author: Marc De Graef
!
!  Description: draw a powder pattern (kinematical)
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!   4/ 3/01 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine DumpPP(xo,yo,np,laL,icnt)

use local
use postscript
use crystalvars

IMPLICIT NONE 

! nref = max number of rings
logical                     :: np
integer(kind=irg),parameter :: nref = 500
integer(kind=irg)           :: i,j,icnt
real(kind=sgl)              :: xo,yo,sc,laL,gmax,leng(nref),PX,PY,locv(nref),grad,w,Vmax
real(kind=sgl),parameter    :: le=3.25,he=2.9375,thr=1.0E-4

! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
 if (np) then
  call PS_newpage(.FALSE.,'Kinematical Zone Axis Patterns')
  call PS_text(5.25,-0.05,'scale bar in reciprocal nm')
  gmax = laL
  call PS_textvar(5.25,PS % psfigheight+0.02,'Camera Constant [nm mm]',gmax)
  call PS_setfont(PSfonts(2),0.15)
  call PS_text(-0.25,PS % psfigheight+0.02,'Structure File : '//cell % fname)
 end if
! draw frame and related stuff
 call PS_setlinewidth(0.012)
 call PS_balloon(xo,yo,le,he,0.0312)
 call PS_setlinewidth(0.001)
 PX = xo+1.8125
 PY = yo+1.0+15.0/32.0
 call PS_circle(PX,PY,1.375)
! frame title
 call PS_setfont(PSfonts(2),0.12)
 call PS_text(xo+0.05,yo+he-0.15,'Powder Pattern')
! scale bar (sc is the conversion factor from nm-1 to inches)
 sc = laL/25.4
 call PS_setlinewidth(0.020)
 call PS_line(xo+0.05,yo+0.06,xo+0.05+5.0*sc,yo+0.06)
 call PS_setfont(PSfonts(2),0.15)
 call PS_text(xo+0.05+2.5*sc,yo+0.10,'5 ')
! scale all reflection intensities by the multiplicity
 leng(1:nref)=0.0
 Vmax = 0.0
 do i=1,icnt-1
  leng(i)=rg(i)
  locv(i)=Vgsave(i)*rnumfam(i)
  if (locv(i).gt.Vmax) Vmax=locv(i)
 end do
! plot origin of reciprocal space 
 call PS_filledcircle(PX,PY,0.05,0.0)
! then plot the diffraction circles 
 do i=1,icnt
  j=icnt+1-i
! get the circle radius and intensity
  grad = leng(j)*sc
  w = locv(j)*0.03/Vmax
! draw circle if radius fits in drawing frame and intensity large enough
  if ((w.gt.0.0001).AND.(grad.le.1.375)) then 
   call PS_setlinewidth(w)
   call PS_circle(PX,PY,grad)
  end if
 end do
end subroutine
! ###################################################################
!
!  subroutine studylist
!
!  Author: Marc De Graef
!
!  Description: determine which zone axis pattern to draw from user input
!               the parameter ppat is returned as true if the 
!               user also requested a powder pattern
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine studylist(list,slect,np,ppat)

use local

IMPLICIT NONE

character(1)                :: list(256)
integer(kind=irg)           :: slect(256),comma(100),hyphen(100),ccnt,hcnt,np,i,j,k,ip,icnt,nd,n,istart,istop
integer(kind=irg),parameter :: nmb(48:57)=(/0,1,2,3,4,5,6,7,8,9/)
logical                     :: ppat

 ccnt = 0
 hcnt = 0
 ppat = .FALSE.
 slect = 0
 comma = 0
 hyphen= 0
 j = 0
! count characters and search for , and -
 do i=1,256
  if (list(i)(1:1).ne.' ') j=j+1
  if (list(i)(1:1).eq.',') then 
   ccnt = ccnt+1
   comma(ccnt)=i
  end if
  if (list(i)(1:1).eq.'-') then 
   hcnt = hcnt+1
   hyphen(hcnt)=i
  end if
 end do 
 ccnt = ccnt+1
 comma(ccnt) = j+1
! interpret the string
 j = 1
 ip = 1
 icnt = 0
 do i=1,ccnt
! is it a range ?
  if (((hyphen(j).lt.comma(i)).and.(hcnt.gt.0)).and.(j.le.hcnt)) then
! yes, it is;  get the first number
   nd = hyphen(j)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   istart = n
   ip = hyphen(j)+1
! and then the second number
   nd = comma(i)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   istop = n
! and fill in the entire range
   do k=istart,istop
    icnt=icnt+1
    slect(icnt)=k
    if (k.eq.0) then
     ppat = .TRUE.
    end if
   end do
   ip = comma(i)+1
   j=j+1
  else
! no, it is not; determine number of digits
   nd = comma(i)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   icnt=icnt+1
   slect(icnt)=n
   if (n.eq.0) then
    ppat = .TRUE.
   end if
   ip = comma(i)+1
  end if
 end do 
 np = icnt
end subroutine
! ###################################################################
!
!  subroutine BWsolve
!
!  Author: Marc De Graef
!
!  Description: Solve the Bloch wave eigenvalue equation for the 
!               N-beam case, using the CGEEV LAPACK 3.0 routine
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!   4/20/01 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine BWsolve(M,W,CG,CGinv,nn,IPIV)

use local
use error

IMPLICIT NONE

integer(kind=irg)    :: INFO, LDA, LDVR, LDVL, LWORK, nn, IPIV(nn), JPIV(nn), i, j, MILWORK
complex(kind=dbl)    :: q,M(nn,nn),VL(nn,nn), CG(nn,nn), W(nn), WORK(4*nn), &
                        MIWORK(nn), CGinv(nn,nn), cW(nn)
real(kind=dbl)       :: RWORK(2*nn)
real(kind=sgl),allocatable :: rW(:)
character            :: JOBVL, JOBVR

intent(IN)  :: nn,M
intent(OUT) :: W,CG,CGinv,IPIV

! first initialize the parameters for the LAPACK CGEEV, CGETRF,
! and CGETRI routines
 JOBVL = 'N'   ! do not compute the left eigenvectors
 JOBVR = 'V'   ! compute the right eigenvectors
 LDA = nn
 LWORK = 4*nn
 LDVL = nn
 LDVR = nn
 MILWORK = nn
! then call the eigenvalue solver
 call ZGEEV(JOBVL,JOBVR,nn,M,LDA,W,VL,LDVL,CG,LDVR,WORK,LWORK,RWORK,INFO)
 if (INFO.ne.0) call FatalError('Error in BWsolve: ','CGEEV return not zero')
! rank the eigenvalues from largest (most positive) to smallest
! (most negative) (use the real part for the ranking); 
! do i=1,nn
!   write (15,*) W(i)
! end do
! write (15,*) '----'
! do i=1,nn
!  do j=1,nn
!   write (15,*) CG(i,j)
!  end do
! write (15,*) '--'
!end do
 cW = W
 allocate(rW(nn))
 rW = real(W)
 IPIV = 0
 JPIV = 0
 call SPSORT(rW,nn,IPIV,-1,INFO)
 do i=1,nn
  W(i) = cW(IPIV(i))
 end do
 CGinv = CG
 do i=1,nn  ! row index
  do j=1,nn ! column index
   CG(i,j) = CGinv(i,IPIV(j))
  end do
 end do 
! make a new copy of CG for the inversion routines
 CGinv = CG
! invert CGinv to get the Bloch wave excitation amplitudes 
 call ZGETRF(nn,nn,CGinv,LDA,JPIV,INFO)
 if (INFO.ne.0) call FatalError('Error in BWsolve: ','CGETRF return not zero')
!write (*,*) 'CGETRF INFO = ',INFO
 call ZGETRI(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
 if (INFO.ne.0) call FatalError('Error in BWsolve: ','CGETRI return not zero')
!write (*,*) 'CGETRI INFO = ',INFO
end subroutine
!
!
!
subroutine BWsolve_test(M,W,CG,CGinv,nn,IPIV)

use local
use error

IMPLICIT NONE

integer    :: ILO, IHI, INFO, LDA, LDVR, LDVL, LWORK, nn, IPIV(nn), JPIV(nn), i, j, MILWORK, IERR
complex    :: q,M(nn,nn),VL(nn,nn), CG(nn,nn), W(nn), WORK(2*(nn**2+2*nn)), &
              MIWORK(nn), CGinv(nn,nn), cW(nn)
real       :: RWORK(2*nn),SCLE(nn),ABNRM,RCONDE(nn),RCONDV(nn),e,slamch
real,allocatable :: rW(:),WR(:),WI(:),CCR(:,:),CCI(:,:)
character  :: JOBVL, JOBVR, BALANC, SENSE

intent(IN)  :: nn,M
intent(OUT) :: W,CG,CGinv,IPIV

! first initialize the parameters for the LAPACK CGEEV, CGETRF,
! and CGETRI routines
 BALANC = 'B'  ! balance and diagonal scaling
 JOBVL = 'V'   ! compute the left eigenvectors
 JOBVR = 'V'   ! compute the right eigenvectors
 SENSE = 'B'   ! reciprocal condition numbers
 LDA = nn
 LWORK = 2*(nn**2+2*nn)
 LDVL = nn
 LDVR = nn
 MILWORK = nn
! then call the eispack eigenvalue solver
! call CGEEVX(BALANC,JOBVL,JOBVR,SENSE,nn,M,LDA,W,VL,LDVL,CG,LDVR, &
!             ILO,IHI,SCLE,ABNRM,RCONDE,RCONDV,WORK,LWORK,RWORK,INFO)
 allocate(WR(nn),WI(nn),CCR(nn,nn),CCI(nn,nn)) 
 CALL EISPACK(nn,nn,real(M),aimag(M),WR,WI,CCR,CCI,IERR)
 W = cmplx(WR,WI)
 CG = cmplx(CCR,CCI)

!if (INFO.ne.0) call FatalError('Error in BWsolve: ','CGEEV return not zero')
  do i=1,nn
    write (15,*) W(i)
  end do
  write (15,*) '----'
  do i=1,nn
   do j=1,nn
    write (15,*) CG(i,j)
   end do
  write (15,*) '--'
 end do

!e = slamch('E')
!write (*,*) 'SLAMCH(E) = ',e
!do i=1,nn
! write (*,*) i,real(W(i)),e*ABNRM/RCONDE(i),e*ABNRM/RCONDV(i)
!end do

! rank the eigenvalues from largest (most positive) to smallest
! (most negative) (use the real part for the ranking); 
 cW = W
 allocate(rW(nn))
 rW = real(W)
 IPIV = 0
 JPIV = 0
 call SPSORT(rW,nn,IPIV,-1,INFO)
 do i=1,nn
  W(i) = cW(IPIV(i))
 end do
 CGinv = CG
 do i=1,nn  ! row index
  do j=1,nn ! column index
   CG(i,j) = CGinv(i,IPIV(j))
  end do
 end do 
! make a new copy of CG for the inversion routines
 CGinv = CG
! invert CGinv to get the Bloch wave excitation amplitudes 
 call CGETRF(nn,nn,CGinv,LDA,JPIV,INFO)
 if (INFO.ne.0) call FatalError('Error in BWsolve: ','CGETRF return not zero')
!write (*,*) 'CGETRF INFO = ',INFO
 call CGETRI(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
 if (INFO.ne.0) call FatalError('Error in BWsolve: ','CGETRI return not zero')
!write (*,*) 'CGETRI INFO = ',INFO
end subroutine
!
! ###################################################################
! 
!  subroutine RankReflections
!
!                                    created: 5/24/01
!  Author: Marc De Graef
!  
!  Description: read a string from standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   5/24/01 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine RankReflections(k,ga,gb,fcnt,srza,iorder)

use local
use multibeams
use symmetryvars
use symmetry
use io
use crystal
use error
use dynamical

IMPLICIT NONE

character(*),parameter :: var1 = 'RankReflections: unable to allocate memory for array '
character(2)        :: srza            ! = SR for systematic row and ZA for zone axis orientation
integer(kind=irg)   :: i,j,ii,jj, &    ! various counters
                       k(3),      &    ! incident beam direction
                       ga(3), gb(3), gn(3), & ! reciprocal lattice vectors
                       il(48), gg(3), iorder,  rfcnt,inm,imm,num,jcnt,iproj,jproj,fcnt,ier
real(kind=sgl)      :: kk(3),Vmod,Vphase,Vpmod,Vpphase,ggg,sga,sgb, &
                       M(2,2),X(2),D ! variables used to decompose reciprocal lattice vectors
logical             :: a,sgg 
logical,allocatable :: z(:,:),zz(:)

intent(IN)          :: k,ga,gb,srza
intent(OUT)         :: fcnt,iorder

 kk = float(k)
 iorder = SG % SYM_numpt
 if (srza.eq.'ZA') then 
  mess = ' The program will list all independent families of reflections of the zone'; call Message("(/A)")
  mess = ' that can be written as linear combinations of the two independent reflections'; call Message("(/A)")
  mess = ' Enter the maximum multiple of ga and gb that should be included [integers] : '; call GetInt(2)
  imm = io_int(1)
  inm = io_int(2)
 else
  mess = ' The program will list all independent families of reflections of the systematic row.'; call Message("(/A)")
  mess = ' Enter the maximum multiple of ga that should be included [I] : '; call GetInt(1)
  inm = io_int(1)
 end if
 fcnt = 0
! systematic row or zone axis ?
 select case (srza)
 case('ZA');
! use a logical array z to exclude family members once one
! member of the family has been evaluated
  if (.not.allocated(z))  allocate(z(-inm:inm,-imm:imm),stat=ier)
  if (ier.ne.0) call FatalError(var1,'z')
  z(-inm:inm,-imm:imm) = .FALSE.
  do i=-inm,inm
   do j=-imm,imm
    if (.not.z(i,j)) then
     fcnt = fcnt + 1
     if (fcnt.gt.numr) call ErrorMess(6)
     gn = i*ga+j*gb
     call CalcFamily(gn,num,'r')
     call GetOrder(kk,il,num,jcnt)
     numfam(fcnt) = jcnt
    do jj = 1,jcnt
      do ii = 1,3
       gg(ii) = itmp(il(jj),ii)
      end do
      if (jj.eq.1) then
       glen(fcnt) = CalcLength(float(gg),'r')
      end if
! determine the components of all family members with
! respect to ga and gb      
      M(1,1) = CalcDot(float(gb),float(gb),'c')
      M(1,2) = -CalcDot(float(ga),float(gb),'c')
      M(2,1) = M(1,2)
      M(2,2) = CalcDot(float(ga),float(ga),'c')
      D = M(1,1)*M(2,2) - M(1,2)*M(2,1)
      X(1) = CalcDot(float(gg),float(ga),'c')
      X(2) = CalcDot(float(gg),float(gb),'c')
      X = matmul(M,X)/D
      iproj = int(X(1))
      jproj = int(X(2))
      if (((iproj.ge.-inm).and.(iproj.le.inm)).and.((jproj.ge.-imm).and.(jproj.le.imm))) then
       z(iproj,jproj) = .TRUE.
      endif
      do ii = 1,3
       family(fcnt,jj,ii) = itmp(il(jj),ii)
      end do
     end do
    end if
   end do
  end do
  deallocate(z)
 case('SR');   ! systematic row case
! In this case we should mostly check for centrosymmetry, since 
! the non-centrosymmetric case may give +g and -g reflections with
! different intensities.
! The most straightforward way to check this is to determine 
! whether or not -g belongs to the family {+g}.  If it does, then
! the systematic row will be symmetric in +-g, else it is not.
  call CalcFamily(ga,num,'r')
  sgg = .FALSE.
  do ii=1,num
   gg(1:3) = itmp(ii,1:3)
   if (sum(gg+ga).eq.0) then 
     sgg = .TRUE.
   end if
  end do
  fcnt = fcnt + 1
  numfam(fcnt) = 1
  glen(fcnt) = 0.0
  family(fcnt,1,1:3) = (/0,0,0/)
  if (sgg) then
! -g does belong to the same family; enumerate all families
   iorder = 2
   do i=1,inm
    fcnt = fcnt + 1
    numfam(fcnt) = 2
    glen(fcnt) = CalcLength(float(i*ga),'r')
    family(fcnt,1,1:3) = i*ga(1:3)
    family(fcnt,2,1:3) = -i*ga(1:3)
   end do
  else
! -g does not belong to the same family
   iorder = 1
   do i=-inm,inm
    fcnt = fcnt + 1
    numfam(fcnt) = 1
    glen(fcnt) = CalcLength(float(i*ga),'r')
    family(fcnt,1,1:3) = i*ga(1:3)
   end do
  end if
 end select   ! systematic row or zone axis
! allocate the variables defined in module beams
 if (.not.allocated(gm)) allocate(gm(fcnt),stat=ier)
 if (ier.ne.0) call FatalError(var1,'gm')
 if (.not.allocated(idx)) allocate(idx(fcnt),stat=ier)
 if (ier.ne.0) call FatalError(var1,'idx')
 if (.not.allocated(al)) allocate(al(fcnt),stat=ier)
 if (ier.ne.0) call FatalError(var1,'al')
 if (.not.allocated(V)) allocate(V(fcnt,4),stat=ier)
 if (ier.ne.0) call FatalError(var1,'V')
 gm=0.0
 gm(1:fcnt) = glen(1:fcnt)
! next, rank the families by increasing |g|
 call SPSORT(gm,fcnt,idx,1,ier)
! and print them out in the correct order, removing the ones that are not
! allowed by lattice centering !  Reflections with zero structure factor,
! but not forbidden by lattice centering may give rise to double diffracted
! beams, so they can not be excluded at this point.
 mess = ' List of independent reflections, ranked by increasing |g|'; call Message("(/A/)")
 rfcnt = 0
 do i=1,fcnt
  gn(1:3) = family(idx(i),1,1:3)
  a = IsGAllowed(gn)
  if (a) then
   call CalcUcg(gn)
   rfcnt=rfcnt+numfam(idx(i))
   al(idx(i)) = .TRUE.
   V(idx(i),1) = rlp%Vmod
   V(idx(i),2) = rlp%Vphase
   V(idx(i),3) = rlp%Vpmod
   V(idx(i),4) = rlp%Vpphase
    write (*,"('(',I3,I3,I3,'); |g| = ',F7.3,'; # = ',I2,'; nref = ',I3,'; |Vg| = ',F8.5,' Volt')") &
        (gn(j),j=1,3),gm(idx(i)),numfam(idx(i)),rfcnt,rlp%Vmod
  else
   al(idx(i)) = .FALSE.
  end if
 end do
end subroutine

!
! ###################################################################
! 
!  subroutine SelectReflections
!
!                                    created: 5/24/01
!  Author: Marc De Graef
!  
!  Description: determine a subset of beams to be included in multi-beam
!               simulation
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   5/24/01 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine SelectReflections(fcnt,rfcnt,ccnt)

use local
use multibeams
use crystalvars
use symmetryvars
use symmetry
use io
use dynamical

IMPLICIT NONE

integer(kind=irg)    :: rfcnt, gn(3), fcnt,i,j,igv, ccnt
real(kind=sgl)       :: gmax

intent(IN) :: fcnt
intent(OUT):: rfcnt,ccnt

! ask for maximum value in terms of |g| or in terms of Vmod
 mess = ' The reciprocal lattice vectors contributing to the computation'; call Message("(/A)")
 mess = ' must now be selected.  You can use a maximum |g| value to '; call Message("(A)")
 mess = ' truncate reciprocal space, or you can use a |V| threshold.'; call Message("(A/)")
! is this a non-symmorphic space group ? If so, add a little note
 if (minval(abs(SGsym - cell % SYM_SGnum)).ne.0) then
  mess = ' The list above may show reflections with zero structure factor because'; call Message("(A)")
  mess = ' of non-symmorphic space group symmetry elements.  Those reflections'; call Message("(A)")
  mess = ' can be incorporated in the simulation only with the maximum |g| criterion.'; call Message("(A/)")
 end if
 mess = 'Use maximum |g| as criterion (1) or threshold |V| (2) :'; call GetInt(1); igv = io_int(1)
 if (igv.eq.1) then
  mess = ' You have selected the maximum |g| criterion.'; call Message("(/A)")
  mess = 'Enter the truncation value for |g| (strictly smaller) : '; call GetReal(1); gmax = io_real(1)
  do i=1,fcnt
   if ((gm(idx(i)).gt.gmax).and.(al(idx(i)))) al(idx(i)) = .FALSE.
  end do 
 else
  mess = ' You have selected the threshold |V| criterion.'; call Message("(/A)")
  mess = 'Enter the minimum value for |V| (strictly larger) : '; call GetReal(1); gmax = io_real(1)
  do i=1,fcnt
   if ((V(idx(i),1).lt.gmax).and.(al(idx(i)))) al(idx(i)) = .FALSE.
  end do 
 endif
 rfcnt = 0
 ccnt = 0
 mess = ' The following reflections will be included :'; call Message("(A/)")
 do i=1,fcnt
  if (al(idx(i))) then 
! place all reflections in the linked reflection list
   do j=1,numfam(idx(i))
    call AddReflection(family(idx(i),j,1:3))
    rltail%famnum = ccnt+1
   end do
! and print list of included families
   gn(1:3) = family(idx(i),1,1:3)
   ccnt = ccnt + 1
   rfcnt=rfcnt+numfam(idx(i))
    write (*,"('(',I3,I3,I3,'); |g| = ',F7.3,'; # = ',I2,'; nref = ',I3,'; |Vg| = ',F8.5,' Volt')") &
        (gn(j),j=1,3),gm(idx(i)),numfam(idx(i)),rfcnt,V(idx(i),1)
  end if
 end do
 mess = 'Total number of reflections contributing to computation : '; oi_int(1) = rfcnt; call WriteInt(1,"(I4)")
end subroutine

!
! ###################################################################
! 
!  subroutine CalcFresnelPropagator
!
!                                    created: 4/16/97
!  Author: Marc De Graef
!  
!  Description: compute the Fresnel propagator (for possibly inclined
!               illumination) and store it in a file
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   4/16/97 MDG 1.0 original
!   9/29/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine CalcFresnelPropagator(beam,dimi,dimj,dz,scl,propname)

use local
use constants
use io
use files

IMPLICIT NONE

real(kind=sgl)                  :: beam(3),b,bm(2),dz,fidim,fjdim,prefac,scl
real(kind=sgl),allocatable      :: idimi(:),jdimj(:)
complex(kind=sgl),allocatable   :: fr(:,:)
integer(kind=irg)               :: dimi,dimj,i,ix,iy
character(20)                   :: propname

INTENT(IN) :: beam,dimi,dimj,dz

  fidim = 1.0/float(dimi)
  fjdim = 1.0/float(dimj)
  prefac = scl*cPi*mLambda*dz
  mess = 'Computing Fresnel propagator'; call Message("(A)")
! normalize the incident beam direction and rescale to the wavevector
  b = sqrt(sum(beam**2))
  bm= beam(1:2)/b/mLambda
  oi_real(1:2) = bm(1:2); mess=' Laue center at '; call WriteReal(2,"(F8.4,',',F8.4)")
! allocate variables
  allocate(fr(dimi,dimj))
  allocate(idimi(dimi),jdimj(dimj))
  idimi = float((/ (i, i=0,dimi-1) /))
  jdimj = float((/ (i, i=0,dimj-1) /))
!
  where(idimi.ge.dimi/2) idimi = idimi-float(dimi)
  where(jdimj.ge.dimj/2) jdimj = jdimj-float(dimj)
!
  idimi = prefac*idimi*fidim
  jdimj = prefac*jdimj*fjdim
!
  idimi = idimi*(idimi + 2.0*bm(1))
  jdimj = jdimj*(jdimj + 2.0*bm(2))
! loop over y axis  
  do iy=1,dimj
! loop over x-axis
    do ix=1,dimi
      fr(ix,iy)=cmplx(cos(idimi(ix)+jdimj(iy)),sin(idimi(ix)+jdimj(iy)))
    end do
  end do
  deallocate(idimi,jdimj)
! and store it in a file
  call SafeOpenFile('d1','UNFORMATTED',propname)
  write (dataunit) dimi,dimj
  write (dataunit) fr
  call SafeCloseFile('d1','KEEP',propname)
  deallocate(fr)
end subroutine

end module

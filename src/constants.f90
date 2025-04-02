!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:constants.f90                                                        !
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
!  FILE: "constants.f90"
!                                    created: 1/5/99 {11:26:50 AM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: physical constants used by various programs; periodic
!  table information; atomic radii; 
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  1/5/99   MDG 1.0 original
!  5/18/01  MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
!
!
!
module constants

use local

IMPLICIT NONE

! various physical constants
! cPi          = pi [dimensionless]
! cLight       = velocity of light [m/s]
! cPlanck      = Planck''s constant [Js]
! cBoltzmann   = Boltmann constant [J/K]
! cPermea      = permeability of vacuum [dimensionless]
! cPermit      = permittivity of vacuum [F/m]
! cCharge      = electron charge [C]
! cRestmass    = electron rest mass [kg]
! cMoment      = electron magnetic moment [J/T]
! cJ2eV        = Joules per eV
real(kind=dbl), parameter :: cPi=3.141592653589D0, cLight = 299792458.D0, &
                             cPlanck = 6.626075D-34, cBoltzmann = 1.380658D-23,  &
                             cPermea = 1.2566371D-6, cPermit = 8.854187817D-12, &
                             cCharge = 1.602177D-19, cRestmass = 9.109389D-31, &
                             cMoment = 9.284770D-24, cJ2eV = 1.602177D-19 

! element symbols
character(2), parameter :: ATOM_sym(98)=(/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
                                          'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca', &
                                          'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                                          'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr', &
                                          'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                                          'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                                          'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                                          'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', &
                                          'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                                          'Pa',' U','Np','Pu','Am','Cm','Bk','Cf'/)

! Shannon-Prewitt ionic radii in nanometer
real(kind=sgl), parameter :: ATOM_SPradii(98)=(/0.010,0.010,0.680,0.300,0.160,0.150,0.148,0.146,0.133,0.500, &
                                                0.098,0.065,0.450,0.380,0.340,0.190,0.181,0.500,0.133,0.940, &
                                                0.068,0.060,0.740,0.690,0.670,0.640,0.630,0.620,0.720,0.740, &
                                                0.113,0.073,0.580,0.202,0.196,0.500,0.148,0.110,0.880,0.770, &
                                                0.067,0.068,0.500,0.500,0.500,0.860,0.126,0.970,0.132,0.930, &
                                                0.076,0.222,0.219,0.500,0.167,0.129,0.104,0.111,0.500,0.108, &
                                                0.050,0.104,0.500,0.970,0.500,0.990,0.500,0.960,0.500,0.940, &
                                                0.050,0.050,0.680,0.600,0.520,0.500,0.500,0.500,0.137,0.112, &
                                                0.140,0.132,0.740,0.230,0.227,0.500,0.175,0.137,0.111,0.990, &
                                                0.090,0.083,0.500,0.108,0.500,0.500,0.500,0.500/)

! atomic (metallic) radii in nanometer (0.100 if not known/applicable)
real(kind=sgl), parameter :: ATOM_MTradii(98)=(/0.100,0.100,0.156,0.112,0.100,0.100,0.100,0.100,0.100,0.100, &
                                                0.191,0.160,0.142,0.100,0.100,0.100,0.100,0.100,0.238,0.196, &
                                                0.160,0.146,0.135,0.128,0.136,0.127,0.125,0.124,0.128,0.137, &
                                                0.135,0.139,0.125,0.116,0.100,0.100,0.253,0.215,0.181,0.160, &
                                                0.147,0.140,0.135,0.133,0.134,0.137,0.144,0.152,0.167,0.158, &
                                                0.161,0.143,0.100,0.100,0.270,0.224,0.187,0.182,0.182,0.181, &
                                                0.100,0.100,0.204,0.178,0.177,0.175,0.176,0.173,0.174,0.193, &
                                                0.173,0.158,0.147,0.141,0.137,0.135,0.135,0.138,0.144,0.155, &
                                                0.171,0.174,0.182,0.168,0.100,0.100,0.100,0.100,0.100,0.180, &
                                                0.163,0.154,0.150,0.164,0.100,0.100,0.100,0.100/)

! atom colors for PostScript drawings
character(3), parameter :: ATOM_color(98)=(/'blu','grn','blu','blu','red','bro','blu','red','grn','grn', &
                                            'blu','pnk','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','grn'/)



end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:                                                                     !
! Copyright (c) 1993, 2001, 2002  A. Weickenmeier/Marc De Graef/CMU             !
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
!  FILE: "others.f90"
!                                    created: 10/13/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: routines written by others but included with written permission.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
module others

contains
! ###################################################################
!
!  subroutine fscatt
!
!  Author: A. Weickenmeier
!
!  Description: Computes the complex scattering amplitude
!
!  History
!
!  converted to Fortran-90 by M. De Graef [5/22/01]
!  includes use of other modules of the CTEMsoft package,
!  in particular io and local
! ###################################################################
complex FUNCTION FSCATT (G,UL,Z,SYMBOL,ACCVLT,ABSFLG,ACCFLG,DWFLG)

! modified: 08.07.93  ! more precise fit to doyle turner table
! modified: 04.03.93  ! expansion in phonon abs. for small g
! 
! Author:
! A.Weickenmeier
! Technische Hochschule Darmstadt
! Hochschulstr. 6
! D-6100 Darmstadt
! Germany 
! bitnet: XLTODA6L@DDATHD21.BITNET
! Tel.: 06151/16-3381
!
! CALCULATES THE COMPLEX SCATTERING AMPLITUDE 
! 
! 
! INPUT:
! ======
!
! G     :   SCATTERING VECTOR G/(4*PI) = S = SIN(THETA)/LAMBDA
! UL    : THE RMS THERMAL DISPLACEMENT OF THE ATOM 
! Z     : ATOM
! ACCVLT: ACCELERATION VOLTAGE IN KV
!   NOT REQUIRED IF ACCFLG = .FALSE.
! ABSFLG: =0 NO ABSORPTION
!         =1 ONLY PHONON
!         =2 ONLY CORE
!         =3 PHONON + CORE
! ACCFLG: IF TRUE THEN F DEPENDS ON THE ACCELERATION VOLTAGE
!         IF FALSE REAL(F) MUST BE MULTIPLIED BY GAMMA,
!         IMAG(F) BY GAMMA^2/K0
!  NOTICE: CORE-CONTRIBUTION TO ABSORPTION DEPENDS ON 
!          ACCELERATION VOLTAGE
!          THEREFORE IF ABSFLG >=1 
!          THEN ACCFLG SHOULD BE .TRUE.
! DWFLG : IF TRUE REAL(F) WILL BE MULTIPLIED WITH THE 
!         DEBYE-WALLER FACTOR
!
! OUTPUT:
! =======
!
! FSCATT: COMPLEX SCATTERING AMPLITUDE
! SYMBOL: SYMBOL OF THE ELEMENT Z
!
! UNITS:
! ======
!
! UNIT OF LENGTH IS ANGSTROEM
! 
! NOTE:
! =====
! 
! THE SCATTERING AMPLITUDES ARE MULTIPLIED BY AN ADDITIONAL 
! FACTOR 4PI.
! 
!**************************************************************************

use local
use io

real(kind=sgl),parameter     :: FOURPI=12.56637062
real(kind=sgl)               :: K0,A(4),B(4)
integer(kind=irg)            :: ABSFLG,Z
CHARACTER(2)                 :: SYMBOL
logical                      :: ACCFLG,DWFLG

! CHECK INPUT
 IF (Z .LT. 1 .OR. Z .GT. 98) THEN
  mess = '  Z = '; oi_int(1) = Z; call WriteInt(1,"(I3,' THIS VALUE IS NOT SUPPLIED!')")
  STOP
 END IF
 IF (UL .LT. 0.) THEN
  mess = '  UL = '; oi_real(1) = UL; call WriteReal(1,"(E13.4,' THIS VALUE IS NOT SUPPLIED!')")
  STOP
 END IF
 IF (G .LT. 0.) THEN
  mess = '  G = '; oi_real(1) = G; call WriteReal(1,"(E13.4,' THIS VALUE IS NOT SUPPLIED!')")
  STOP
 END IF
 IF (ACCVLT .LT. 0.) THEN
  mess = '  ACCVLT = '; oi_real(1) = ACCVLT; call WriteReal(1,"(E13.4,' THIS VALUE IS NOT SUPPLIED!')")
  STOP
 END IF

 IF (DWFLG) THEN
!    CALCULATE DEBYE-WALLER FACTOR
  DEWA  = EXP(-.5*UL*UL*G*G)
 ELSE
  DEWA = 1.
 END IF

! GET FITTING COEFFICIENTS
 CALL GETWK (Z,SYMBOL,A,B)

 S      = G / FOURPI
 FREAL  = FOURPI * DEWA * WEKO (A,B,S)

 IF (ABSFLG .EQ. 0) THEN
  FIMA   = 0.
 ELSE IF (ABSFLG .EQ. 1) THEN
       FIMA   = FPHON (G,UL,A,B) 
      ELSE IF (ABSFLG .EQ. 2) THEN
            FIMA   = FCORE (G,Z,ACCVLT) * DEWA
            ELSE IF (ABSFLG .EQ. 3) THEN
                  FIMA   = FCORE (G,Z,ACCVLT) * DEWA + FPHON (G,UL,A,B) 
 END IF
 IF (ACCFLG.EQV..TRUE.) THEN
!    CALCULATE WAVENUMBER AND GAMMA
  K0   = .5068 * SQRT(1022.*ACCVLT + ACCVLT*ACCVLT)
  GAMMA = (ACCVLT+511.) / 511.
  FREAL = FREAL * GAMMA
  FIMA  = FIMA  * GAMMA * GAMMA / K0
 END IF

 FSCATT = CMPLX( FREAL,FIMA )

END FUNCTION
!*******************************************************************
FUNCTION WEKO (A,B,S)

! ELECTRON SCATTERING AMPLITUDE F(S)

! UPDATE 07.07.93
! NOW 4 A AND 4 B

use local

REAL(kind=sgl)   :: A(4), B(4)

 WEKO=0.
 IF (S .GT. 0.) THEN
    S2 = 1./(S*S)
 END IF

 DO I=1,4
  ARGU = B(I)*S*S
  IF (ARGU .LT. .1) THEN
   WEKO = WEKO + A(I)*B(I) * (1.-.5*ARGU)
  ELSE IF (ARGU .GT. 20.) THEN
        WEKO = WEKO + A(I)*S2
       ELSE
        WEKO = WEKO + A(I)*(1.-EXP(-ARGU))*S2
  END IF
 END DO
END FUNCTION

! ***************************************************************
real FUNCTION FPHON (G,UL,A,B)

use local

real(kind=sgl),parameter  :: FOURPI = 12.56636, FP2= FOURPI*FOURPI
real(kind=sgl)            :: A(4),B(4),A1(4),B1(4)

 U2 = UL*UL

 DO I=1,4
  A1(I) = A(I) * FP2 
 END DO
 DO I=1,4
  B1(I) = B(I) / FP2
 END DO
 
 FPHON = 0.
 G2    = G*G
 DEWA  = EXP(-.5*U2*G2)

 DO J=1,4
  FPHON = FPHON + A1(J)*A1(J)*(DEWA * RI1(B1(J),B1(J),G) - RI2(B1(J),B1(J),G,UL))
  DO I=1,J-1
   FPHON = FPHON + 2.*A1(J)*A1(I)*(DEWA * RI1(B1(I),B1(J),G)-RI2(B1(I),B1(J),G,UL))
  END DO
 END DO
END FUNCTION

!*****************************************************************
real FUNCTION RI1 (BI,BJ,G)

! ERSTES INTEGRAL FUER DIE ABSORPTIONSPOTENTIALE

use local

real(kind=sgl),parameter  :: PI=3.1415927, C =0.5772157


 G2   = G*G

! IST DIE ASYMPTOTISCHE ENTWICKLUNG ANWENDBAR?

 EPS = MAX( BI,BJ )
 EPS = EPS * G2

 IF (EPS .LE. 0.1) THEN
  RI1 = BI * LOG( (BI+BJ)/BI ) + BJ * LOG( (BI+BJ)/BJ )
  RI1 = RI1 * PI
  IF (G .EQ. 0.) RETURN
  BI2 = BI*BI
  BJ2 = BJ*BJ
  TEMP = .5*BI2*LOG( BI/(BI+BJ) ) + .5*BJ2*LOG( BJ/(BI+BJ) )
  TEMP = TEMP + .75 * (BI2+BJ2) - .25 * (BI+BJ)*(BI+BJ)
  TEMP = TEMP - .5 * (BI-BJ)*(BI-BJ)
  RI1  = RI1 + PI * G2 * TEMP
  RETURN
 END IF

 BIG2 = BI*G2
 BJG2 = BJ*G2
 RI1 = 2.*C + LOG(BIG2) + LOG(BJG2) - 2.*EI(-BI*BJ*G2/(BI+BJ))

 X1  = BIG2
 X2  = BIG2*BI/(BI+BJ)
 X3  = BIG2
 RI1 = RI1 + RIH1(X1,X2,X3)

 X1  = BJG2
 X2  = BJG2*BJ/(BI+BJ)
 X3  = BJG2
 RI1 = RI1 + RIH1(X1,X2,X3)
 RI1 = RI1 * PI / G2

END FUNCTION

!**************************************************************************

real FUNCTION RI2 (BI,BJ,G,U)

! ZWEITES INTEGRAL FUER DIE ABSORPTIONSPOTENTIALE

use local

real(kind=sgl),parameter  :: PI=3.1415927

 U2 = U*U
 U22 = .5*U2
 G2   = G*G
 BIUH = BI + .5*U2
 BJUH = BJ + .5*U2
 BIU  = BI + U2
 BJU  = BJ + U2

! IST DIE ASYMPTOTISCHE ENTWICKLUNG ANWENDBAR?
 EPS = MAX( BI,BJ,U2 )
 EPS = EPS * G2

 IF (EPS .LE. 0.1) THEN
  RI2 = (BI+U2) * LOG( (BI+BJ+U2)/(BI+U2) )
  RI2 = RI2 + BJ * LOG( (BI+BJ+U2)/(BJ+U2) )
  RI2 = RI2 + U2 * LOG( U2/(BJ+U2) )
  RI2 = RI2 * PI
  IF (G .EQ. 0.) RETURN
  TEMP = .5 * U22 * U22 * LOG( BIU*BJU/(U2*U2) )
  TEMP = TEMP + .5 * BIUH * BIUH * LOG( BIU/(BIUH+BJUH) ) 
  TEMP = TEMP + .5 * BJUH * BJUH * LOG( BJU/(BIUH+BJUH) ) 
  TEMP = TEMP + .25 * BIU * BIU + .5 * BI * BI
  TEMP = TEMP + .25 * BJU * BJU + .5 * BJ * BJ
  TEMP = TEMP - .25 * (BIUH+BJUH) * (BIUH+BJUH)
  TEMP = TEMP - .5 * ( (BI*BIU-BJ*BJU)/(BIUH+BJUH) )**2
  TEMP = TEMP - U22*U22
  RI2  = RI2 + PI * G2 * TEMP
  RETURN
 END IF

 RI2 = EI( -.5*U2*G2*BIUH/BIU ) + EI( -.5*U2*G2*BJUH/BJU )
 RI2 = RI2 - EI( -BIUH*BJUH*G2/(BIUH+BJUH) ) - EI( -.25*U2*G2 )
 RI2 = 2.*RI2
 X1  = .5*U2*G2
 X2  = .25*U2*G2
 X3  = .25*U2*U2*G2/BIU
 RI2 = RI2 + RIH1(X1,X2,X3)

 X1  = .5*U2*G2
 X2  = .25*U2*G2
 X3  = .25*U2*U2*G2/BJU
 RI2 = RI2 + RIH1(X1,X2,X3)

 X1  = BIUH*G2
 X2  = BIUH*BIUH*G2/(BIUH+BJUH)
 X3  = BIUH*BIUH*G2/BIU
 RI2 = RI2 + RIH1(X1,X2,X3)

 X1  = BJUH*G2
 X2  = BJUH*BJUH*G2/(BIUH+BJUH)
 X3  = BJUH*BJUH*G2/BJU
 RI2 = RI2 + RIH1(X1,X2,X3)

 RI2 = RI2 * PI / G2

END FUNCTION

!**************************************************************************
real  FUNCTION RIH1 (X1,X2,X3)

! WERTET DEN AUSDRUCK EXP(-X1) * ( EI(X2)-EI(X3) ) AUS

 IF (X2 .LE. 20.  .AND.  X3 .LE. 20.) THEN
  RIH1 = EXP(-X1) * ( EI(X2)-EI(X3) )
  RETURN
 END IF

 IF (X2 .GT. 20) THEN
  RIH1 = EXP(X2-X1)*RIH2(X2)/X2
 ELSE 
  RIH1 = EXP(-X1)*EI(X2)
 END IF

 IF (X3 .GT. 20) THEN
  RIH1 = RIH1 - EXP(X3-X1)*RIH2(X3)/X3
 ELSE 
  RIH1 = RIH1 - EXP(-X1)*EI(X3)
 END IF
END FUNCTION

!**************************************************************************
real FUNCTION RIH2 (X)

! WERTET X*EXP(-X)*EI(X) AUS FUER GROSSE X
! DURCH INTERPOLATION DER TABELLE ... AUS ABRAMOWITZ

use local

real(kind=sgl),parameter :: F(0:20) =(/1.000000,1.005051,1.010206,1.015472,1.020852, &
                                       1.026355,1.031985,1.037751,1.043662,1.049726, &
                                       1.055956,1.062364,1.068965,1.075780,1.082830, &
                                       1.090140,1.097737,1.105647,1.113894,1.122497, &
                                       1.131470/)

 X1 = 1./X
 I  = INT( 200.*X1 )
 I1 = I+1

 RIH2 = F(I) + 200.*( F(I1)-F(I) ) * ( X1-.5E-3*REAL(I) )
END FUNCTION

!**************************************************************************
real FUNCTION RIH3 (X)

! WERTET DEN AUSDRUCK EXP(-X) * EI(X) AUS

 IF (X .LE. 20.) THEN
  RIH3 = EXP(-X) * EI(X)
  RETURN
 ELSE
  RIH3 = RIH2(X)/X
 END IF
END FUNCTION

! ***********************************************************
real FUNCTION EI (X)

! EXPONENTIALINTEGRAL
!   GETESTET -60 < X < 60

use local
use io

real(kind=sgl),parameter :: A1=8.57332,A2=18.05901,A3=8.63476,A4=.26777, &
                            B1=9.57332,B2=25.63295,B3=21.09965,B4=3.95849

 IF (X .GT. 60.) THEN 
  mess = '>>> EI NICHT GETESTET FUER X= '; oi_real(1) = X; call WriteInt(1,"(f10.5,' <<<<')")
  STOP
 END IF

 IF (X .LT. -60.) THEN
  EI = 0.
  RETURN
 END IF

 IF (X .LT. -1.) THEN
!    ABRAMOWITZ (5.1.56)
  XP = ABS(X)
  EI = -( A4+XP*(A3+XP*(A2+XP*(A1+XP))) )/( B4+XP*(B3+XP*(B2+XP*(B1+XP))) ) * EXP(-XP)/XP
  RETURN
 END IF
    
 EI   = .577216 + LOG( ABS(X) )

 I    = 1
 SI   = X
 SUMM = SI
 do while (ABS(SI/X) .GT. 1.E-6)
  RI   = REAL(I)
  RII  = RI + 1.
  SI   = SI * X * RI/(RII*RII)
  SUMM = SUMM + SI
  I    = I+1    
 end do
 EI = EI + SUMM
END FUNCTION

!**************************************************************************
SUBROUTINE GETWK (Z,SYMBOL,A,B)

! UPDATE 07.07.93
! NOW 4 A AND 4 B

use local

! DATEN VON UND FUER DIE AUFRUFENDE ROUTINE
integer       :: Z
real(kind=sgl):: A(4),B(4)
character(2)  :: SYMBOL 

! DATEN, DIE NUR INTERN BENOETIGT WERDEN
character(2),parameter :: SY(98) = (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
                                     'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
                                     'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                                     'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
                                     'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                                     'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                                     'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                                     'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
                                     'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                                     'Pa','U ','Np','Pu','Am','Cm','Bk','Cf'/)

REAL(kind=sgl),parameter :: AA(4,98) = reshape((/ 0.00427, 0.00957, 0.00802, 0.00209, &
                                   0.01217, 0.02616,-0.00884, 0.01841, &
                                   0.00251, 0.03576, 0.00988, 0.02370, &
                                   0.01596, 0.02959, 0.04024, 0.01001, &
                                   0.03652, 0.01140, 0.05677, 0.01506, &
                                   0.04102, 0.04911, 0.05296, 0.00061, &
                                   0.04123, 0.05740, 0.06529, 0.00373, &
                                   0.03547, 0.03133, 0.10865, 0.01615, &
                                   0.03957, 0.07225, 0.09581, 0.00792, &
                                   0.02597, 0.02197, 0.13762, 0.05394, &
                                   0.03283, 0.08858, 0.11688, 0.02516, &
                                   0.03833, 0.17124, 0.03649, 0.04134, &
                                   0.04388, 0.17743, 0.05047, 0.03957, &
                                   0.03812, 0.17833, 0.06280, 0.05605, &
                                   0.04166, 0.17817, 0.09479, 0.04463, &
                                   0.04003, 0.18346, 0.12218, 0.03753, &
                                   0.04245, 0.17645, 0.15814, 0.03011, &
                                   0.05011, 0.16667, 0.17074, 0.04358, &
                                   0.04058, 0.17582, 0.20943, 0.02922, &
                                   0.04001, 0.17416, 0.20986, 0.05497, &
                                   0.09685, 0.14777, 0.20981, 0.04852, &
                                   0.06667, 0.17356, 0.22710, 0.05957, &
                                   0.05118, 0.16791, 0.26700, 0.06476, &
                                   0.03204, 0.18460, 0.30764, 0.05052, &
                                   0.03866, 0.17782, 0.31329, 0.06898, &
                                   0.05455, 0.16660, 0.33208, 0.06947, &
                                   0.05942, 0.17472, 0.34423, 0.06828, &
                                   0.06049, 0.16600, 0.37302, 0.07109, &
                                   0.08034, 0.15838, 0.40116, 0.05467, &
                                   0.02948, 0.19200, 0.42222, 0.07480, &
                                   0.16157, 0.32976, 0.18964, 0.06148, &
                                   0.16184, 0.35705, 0.17618, 0.07133, &
                                   0.06190, 0.18452, 0.41600, 0.12793, &
                                   0.15913, 0.41583, 0.13385, 0.10549, &
                                   0.16514, 0.41202, 0.12900, 0.13209, &
                                   0.15798, 0.41181, 0.14254, 0.14987, &
                                   0.16535, 0.44674, 0.24245, 0.03161, &
                                   0.16039, 0.44470, 0.24661, 0.05840, &
                                   0.16619, 0.44376, 0.25613, 0.06797, &
                                   0.16794, 0.44505, 0.27188, 0.07313, &
                                   0.16552, 0.45008, 0.30474, 0.06161, &
                                   0.17327, 0.44679, 0.32441, 0.06143, &
                                   0.16424, 0.45046, 0.33749, 0.07766, &
                                   0.18750, 0.44919, 0.36323, 0.05388, &
                                   0.16081, 0.45211, 0.40343, 0.06140, &
                                   0.16599, 0.43951, 0.41478, 0.08142, &
                                   0.16547, 0.44658, 0.45401, 0.05959, &
                                   0.17154, 0.43689, 0.46392, 0.07725, &
                                   0.15752, 0.44821, 0.48186, 0.08596, &
                                   0.15732, 0.44563, 0.48507, 0.10948, &
                                   0.16971, 0.42742, 0.48779, 0.13653, &
                                   0.14927, 0.43729, 0.49444, 0.16440, &
                                   0.18053, 0.44724, 0.48163, 0.15995, &
                                   0.13141, 0.43855, 0.50035, 0.22299, &
                                   0.31397, 0.55648, 0.39828, 0.04852, &
                                   0.32756, 0.53927, 0.39830, 0.07607, &
                                   0.30887, 0.53804, 0.42265, 0.09559, &
                                   0.28398, 0.53568, 0.46662, 0.10282, &
                                   0.35160, 0.56889, 0.42010, 0.07246, &
                                   0.33810, 0.58035, 0.44442, 0.07413, &
                                   0.35449, 0.59626, 0.43868, 0.07152, &
                                   0.35559, 0.60598, 0.45165, 0.07168, &
                                   0.38379, 0.64088, 0.41710, 0.06708, &
                                   0.40352, 0.64303, 0.40488, 0.08137, &
                                   0.36838, 0.64761, 0.47222, 0.06854, &
                                   0.38514, 0.68422, 0.44359, 0.06775, &
                                   0.37280, 0.67528, 0.47337, 0.08320, &
                                   0.39335, 0.70093, 0.46774, 0.06658, &
                                   0.40587, 0.71223, 0.46598, 0.06847, &
                                   0.39728, 0.73368, 0.47795, 0.06759, &
                                   0.40697, 0.73576, 0.47481, 0.08291, &
                                   0.40122, 0.78861, 0.44658, 0.08799, &
                                   0.41127, 0.76965, 0.46563, 0.10180, &
                                   0.39978, 0.77171, 0.48541, 0.11540, &
                                   0.39130, 0.80752, 0.48702, 0.11041, &
                                   0.40436, 0.80701, 0.48445, 0.12438, &
                                   0.38816, 0.80163, 0.51922, 0.13514, &
                                   0.39551, 0.80409, 0.53365, 0.13485, &
                                   0.40850, 0.83052, 0.53325, 0.11978, &
                                   0.40092, 0.85415, 0.53346, 0.12747, &
                                   0.41872, 0.88168, 0.54551, 0.09404, &
                                   0.43358, 0.88007, 0.52966, 0.12059, &
                                   0.40858, 0.87837, 0.56392, 0.13698, &
                                   0.41637, 0.85094, 0.57749, 0.16700, &
                                   0.38951, 0.83297, 0.60557, 0.20770, &
                                   0.41677, 0.88094, 0.55170, 0.21029, &
                                   0.50089, 1.00860, 0.51420, 0.05996, &
                                   0.47470, 0.99363, 0.54721, 0.09206, &
                                   0.47810, 0.98385, 0.54905, 0.12055, &
                                   0.47903, 0.97455, 0.55883, 0.14309, &
                                   0.48351, 0.98292, 0.58877, 0.12425, &
                                   0.48664, 0.98057, 0.61483, 0.12136, &
                                   0.46078, 0.97139, 0.66506, 0.13012, &
                                   0.49148, 0.98583, 0.67674, 0.09725, &
                                   0.50865, 0.98574, 0.68109, 0.09977, &
                                   0.46259, 0.97882, 0.73056, 0.12723, &
                                   0.46221, 0.95749, 0.76259, 0.14086, &
                                   0.48500, 0.95602, 0.77234, 0.13374/),(/4,98/))
                         

REAL(kind=sgl),parameter :: BB(4,98) = reshape((/ 4.17218, 16.05892, 26.78365, 69.45643, &
                                   1.83008,  7.20225, 16.13585, 18.75551, &
                                   0.02620,  2.00907, 10.80597,130.49226, &
                                   0.38968,  1.99268, 46.86913,108.84167, &
                                   0.50627,  3.68297, 27.90586, 74.98296, &
                                   0.41335, 10.98289, 34.80286,177.19113, &
                                   0.29792,  7.84094, 22.58809, 72.59254, &
                                   0.17964,  2.60856, 11.79972, 38.02912, &
                                   0.16403,  3.96612, 12.43903, 40.05053, &
                                   0.09101,  0.41253,  5.02463, 17.52954, &
                                   0.06008,  2.07182,  7.64444,146.00952, &
                                   0.07424,  2.87177, 18.06729, 97.00854, &
                                   0.09086,  2.53252, 30.43883, 98.26737, &
                                   0.05396,  1.86461, 22.54263, 72.43144, &
                                   0.05564,  1.62500, 24.45354, 64.38264, &
                                   0.05214,  1.40793, 23.35691, 53.59676, &
                                   0.04643,  1.15677, 19.34091, 52.88785, &
                                   0.07991,  1.01436, 15.67109, 39.60819, &
                                   0.03352,  0.82984, 14.13679,200.97722, &
                                   0.02289,  0.71288, 11.18914,135.02390, &
                                   0.12527,  1.34248, 12.43524,131.71112, &
                                   0.05198,  0.86467, 10.59984,103.56776, &
                                   0.03786,  0.57160,  8.30305, 91.78068, &
                                   0.00240,  0.44931,  7.92251, 86.64058, &
                                   0.01836,  0.41203,  6.73736, 76.30466, &
                                   0.03947,  0.43294,  6.26864, 71.29470, &
                                   0.03962,  0.43253,  6.05175, 68.72437, &
                                   0.03558,  0.39976,  5.36660, 62.46894, &
                                   0.05475,  0.45736,  5.38252, 60.43276, &
                                   0.00137,  0.26535,  4.48040, 54.26088, &
                                   0.10455,  2.18391,  9.04125, 75.16958, &
                                   0.09890,  2.06856,  9.89926, 68.13783, &
                                   0.01642,  0.32542,  3.51888, 44.50604, &
                                   0.07669,  1.89297, 11.31554, 46.32082, &
                                   0.08199,  1.76568,  9.87254, 38.10640, &
                                   0.06939,  1.53446,  8.98025, 33.04365, &
                                   0.07044,  1.59236, 17.53592,215.26198, &
                                   0.06199,  1.41265, 14.33812,152.80257, &
                                   0.06364,  1.34205, 13.66551,125.72522, &
                                   0.06565,  1.25292, 13.09355,109.50252, &
                                   0.05921,  1.15624, 13.24924, 98.69958, &
                                   0.06162,  1.11236, 12.76149, 90.92026, &
                                   0.05081,  0.99771, 11.28925, 84.28943, &
                                   0.05120,  1.08672, 12.23172, 85.27316, &
                                   0.04662,  0.85252, 10.51121, 74.53949, &
                                   0.04933,  0.79381,  9.30944, 41.17414, &
                                   0.04481,  0.75608,  9.34354, 67.91975, &
                                   0.04867,  0.71518,  8.40595, 64.24400, &
                                   0.03672,  0.64379,  7.83687, 73.37281, &
                                   0.03308,  0.60931,  7.04977, 64.83582, &
                                   0.04023,  0.58192,  6.29247, 55.57061, &
                                   0.02842,  0.50687,  5.60835, 48.28004, &
                                   0.03830,  0.58340,  6.47550, 47.08820, &
                                   0.02097,  0.41007,  4.52105, 37.18178, &
                                   0.07813,  1.45053, 15.05933,199.48830, &
                                   0.08444,  1.40227, 13.12939,160.56676, &
                                   0.07206,  1.19585, 11.55866,127.31371, &
                                   0.05717,  0.98756,  9.95556,117.31874, &
                                   0.08249,  1.43427, 12.37363,150.55968, &
                                   0.07081,  1.31033, 11.44403,144.17706, &
                                   0.07442,  1.38680, 11.54391,143.72185, &
                                   0.07155,  1.34703, 11.00432,140.09138, &
                                   0.07794,  1.55042, 11.89283,142.79585, &
                                   0.08508,  1.60712, 11.45367,116.64063, &
                                   0.06520,  1.32571, 10.16884,134.69034, &
                                   0.06850,  1.43566, 10.57719,131.88972, &
                                   0.06264,  1.26756,  9.46411,107.50194, &
                                   0.06750,  1.35829,  9.76480,127.40374, &
                                   0.06958,  1.38750,  9.41888,122.10940, &
                                   0.06574,  1.31578,  9.13448,120.98209, &
                                   0.06517,  1.29452,  8.67569,100.34878, &
                                   0.06213,  1.30860,  9.18871, 91.20213, &
                                   0.06292,  1.23499,  8.42904, 77.59815, &
                                   0.05693,  1.15762,  7.83077, 67.14066, &
                                   0.05145,  1.11240,  8.33441, 65.71782, &
                                   0.05573,  1.11159,  8.00221, 57.35021, &
                                   0.04855,  0.99356,  7.38693, 51.75829, &
                                   0.04981,  0.97669,  7.38024, 44.52068, &
                                   0.05151,  1.00803,  8.03707, 45.01758, &
                                   0.04693,  0.98398,  7.83562, 46.51474, &
                                   0.05161,  1.02127,  9.18455, 64.88177, &
                                   0.05154,  1.03252,  8.49678, 58.79463, &
                                   0.04200,  0.90939,  7.71158, 57.79178, &
                                   0.04661,  0.87289,  6.84038, 51.36000, &
                                   0.04168,  0.73697,  5.86112, 43.78613, &
                                   0.04488,  0.83871,  6.44020, 43.51940, &
                                   0.05786,  1.20028, 13.85073,172.15909, &
                                   0.05239,  1.03225, 11.49796,143.12303, &
                                   0.05167,  0.98867, 10.52682,112.18267, &
                                   0.04931,  0.95698,  9.61135, 95.44649, &
                                   0.04748,  0.93369,  9.89867,102.06961, &
                                   0.04660,  0.89912,  9.69785,100.23434, &
                                   0.04323,  0.78798,  8.71624, 92.30811, &
                                   0.04641,  0.85867,  9.51157,111.02754, &
                                   0.04918,  0.87026,  9.41105,104.98576, &
                                   0.03904,  0.72797,  8.00506, 86.41747, &
                                   0.03969,  0.68167,  7.29607, 75.72682, &
                                   0.04291,  0.69956,  7.38554, 77.18528/),(/4,98/))

 SYMBOL = SY(Z)
 DO K=1,4
    A(K) = AA(K,Z)
    B(K) = BB(K,Z)
 END DO
END SUBROUTINE
 
!*******************************************************************
real FUNCTION FCORE (G,Z,ACCVLT)

! modified: 13.01.1992
! 
! Author:
! A.Weickenmeier
! Technische Hochschule Darmstadt
! Hochschulstr. 6
! D-6100 Darmstadt
! Germany 
! bitnet: XLTODA6L@DDATHD21.BITNET
! Tel.: 06151/16-3381
!
! CALCULATES APPROXIMATELY THE CORE PART OF THE APSORTION
! USING THE APROXIMATION FOR THE MIXED DYNAMIC FORM FACTOR
! AND THE X-RAY SCATTERING AMPLITUDE DUE TO
! H. ROSE, OPTIK 45  NO.2  (1976) P. 139-158.
! 
! 
! INPUT:
! ======
!
! G:   SCATTERING VECTOR 4*PI*G = S = SIN(THETA)/LAMBDA
! ACCVLT: ACCELERATION VOLTAGE IN KV
! Z: ATOM
!
! OUTPUT:
! =======
!
! FCORE : CORE PART OF ABSORPTION POTENTIAL
!
! UNITS:
! ======
!
! UNIT OF LENGTH IS ANGSTROEM
! 
!**************************************************************************

use local
use io

real(kind=sgl),parameter  :: PI=3.1415927,A0=.5289
integer                   :: Z
real(kind=sgl)            :: K0,K2,KAPPA

! CHECK INPUT
 IF (Z .LT. 1) THEN
  mess = '  Z = '; oi_int(1) = Z; call WriteInt(1,"(I3,' THIS VALUE IS NOT SUPPLIED!')")
  STOP
 END IF
 IF (G .LT. 0.) THEN
  mess = '  G = '; oi_real(1) = G; call WriteReal(1,"(E13.4,' THIS VALUE IS NOT SUPPLIED!')")
  STOP
 END IF
 IF (ACCVLT .LT. 0.) THEN
  mess = '  ACCVLT = '; oi_real(1) = ACCVLT; call WriteReal(1,"(E13.4,' THIS VALUE IS NOT SUPPLIED!')")
  STOP
 END IF

! CALCULATE WAVENUMBER 
 K0   = .5068 * SQRT(1022.*ACCVLT + ACCVLT*ACCVLT)
! CALCULATE CHARACTERISTIC ENERGY LOSS AND ANGLE
 DE = REAL (6*Z) * 1.E-3
 THETAE = DE/(2.*ACCVLT) * (2.*ACCVLT+1022.)/(ACCVLT+1022.)
! SCREENING PARAMETER OF YUKAWA-POTENTIAL
 R = 0.885 * A0 / REAL(Z)**.333333
! CALCULATE NORMALISING ANGLE
 TA = 1./(K0 * R)
! CALCULATE BRAGG ANGLE
 TB = G / (2.*K0)
! NORMALIZE
 OMEGA = 2.*TB/TA
 KAPPA = THETAE/TA

 K2 = KAPPA*KAPPA
 O2 = OMEGA*OMEGA

 X1 = OMEGA/((1.+O2)*SQRT(O2+4.*K2))*LOG((OMEGA + SQRT( O2+4.*K2 ))/(2.*KAPPA))
 X2 = 1./SQRT((1.+O2)*(1.+O2)+4.*K2*O2)*LOG((1.+2.*K2+O2+SQRT((1.+O2)*(1.+O2)+4.*K2*O2))/(2.*KAPPA*SQRT(1.+K2))) 
 IF (OMEGA .GT. 1E-2) THEN
  X3 = 1./(OMEGA*SQRT( O2+4.*(1.+K2)))* LOG((OMEGA + SQRT( O2+4.*(1.+K2)))/(2.*SQRT(1.+K2))) 
 ELSE
  X3 = 1./(4.*(1.+K2))
 END IF

 HI = REAL(2*Z)/(TA*TA) * (-X1+X2-X3)

 FCORE = 4./(A0*A0) * 2.*PI/(K0*K0) * HI
END FUNCTION
!
! END OF WEICKENMEIER-KOHL SOURCE CODE
!

end module

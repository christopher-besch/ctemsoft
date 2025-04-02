!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:symmetry.f90                                                         !
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
! 
!  FILE: "symmetry.f90"
!                                    created: 1/5/99 {11:26:51 AM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: all symmetry and space group routines
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  1/5/99   MDG 1.0 original
!  5/19/01  MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
!

module symmetry

use local
use symmetryvars

contains

!  Description: This module contains all the space group routines
! 
!    subroutine SYM_fillgen(t,isgn)
!    subroutine MakeGenerators
!    subroutine matrixmult(k1,k2)
!    subroutine GenerateSymmetry
!    subroutine CalcFamily(ind,num,space)
!    subroutine CalcOrbit(m,n,ctmp)
!    subroutine CalcPositions(switch)
!    logical function isitnew(nsym)
!
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
!
! ###################################################################
! 
!  subroutine SYM_fillgen
!
!  Author: Marc De Graef
!  
!  Description: create a generator matrix and translation vector
!               based on the 4-character code strings
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine SYM_fillgen(t,isgn)

use local

IMPLICIT NONE

character(1)           :: t(4)
integer(kind=irg)      :: j,isgn
real(kind=sgl)         :: sgn

intent(IN)             :: t,isgn

! forward or reverse translation ?
 sgn=float(isgn)
! first fill the array with zeroes and a 1 at 4,4
 SG%SYM_c(1:4,1:4) = 0.0_dbl
 SG%SYM_c(4,4) = 1.0_dbl
! then check for the particular matrix
 select case (t(1))
  case ('a'); SG%SYM_c(1,1) = 1.0_dbl; SG%SYM_c(2,2) = 1.0_dbl; SG%SYM_c(3,3) = 1.0_dbl
  case ('b'); SG%SYM_c(1,1) =-1.0_dbl; SG%SYM_c(2,2) =-1.0_dbl; SG%SYM_c(3,3) = 1.0_dbl
  case ('c'); SG%SYM_c(1,1) =-1.0_dbl; SG%SYM_c(2,2) = 1.0_dbl; SG%SYM_c(3,3) =-1.0_dbl
  case ('d'); SG%SYM_c(1,3) = 1.0_dbl; SG%SYM_c(2,1) = 1.0_dbl; SG%SYM_c(3,2) = 1.0_dbl
  case ('e'); SG%SYM_c(1,2) = 1.0_dbl; SG%SYM_c(2,1) = 1.0_dbl; SG%SYM_c(3,3) =-1.0_dbl
  case ('f'); SG%SYM_c(1,2) =-1.0_dbl; SG%SYM_c(2,1) =-1.0_dbl; SG%SYM_c(3,3) =-1.0_dbl
  case ('g'); SG%SYM_c(1,2) =-1.0_dbl; SG%SYM_c(2,1) = 1.0_dbl; SG%SYM_c(3,3) = 1.0_dbl
  case ('h'); SG%SYM_c(1,1) =-1.0_dbl; SG%SYM_c(2,2) =-1.0_dbl; SG%SYM_c(3,3) =-1.0_dbl
  case ('i'); SG%SYM_c(1,1) = 1.0_dbl; SG%SYM_c(2,2) = 1.0_dbl; SG%SYM_c(3,3) =-1.0_dbl
  case ('j'); SG%SYM_c(1,1) = 1.0_dbl; SG%SYM_c(2,2) =-1.0_dbl; SG%SYM_c(3,3) = 1.0_dbl
  case ('k'); SG%SYM_c(1,2) =-1.0_dbl; SG%SYM_c(2,1) =-1.0_dbl; SG%SYM_c(3,3) = 1.0_dbl
  case ('l'); SG%SYM_c(1,2) = 1.0_dbl; SG%SYM_c(2,1) = 1.0_dbl; SG%SYM_c(3,3) = 1.0_dbl
  case ('m'); SG%SYM_c(1,2) = 1.0_dbl; SG%SYM_c(2,1) =-1.0_dbl; SG%SYM_c(3,3) =-1.0_dbl
  case ('n'); SG%SYM_c(1,2) =-1.0_dbl; SG%SYM_c(2,1) = 1.0_dbl; SG%SYM_c(2,2) =-1.0_dbl; SG%SYM_c(3,3) = 1.0_dbl
 end select
! then fill in the translational component
 do j=2,4 
  select case (t(j))
   case('A'); SG%SYM_c(j-1,4) = sgn/6.0_dbl
   case('B'); SG%SYM_c(j-1,4) = sgn/4.0_dbl
   case('C'); SG%SYM_c(j-1,4) = sgn/3.0_dbl
   case('D'); SG%SYM_c(j-1,4) = sgn/2.0_dbl
   case('E'); SG%SYM_c(j-1,4) = sgn*2.0_dbl/3.0_dbl
   case('F'); SG%SYM_c(j-1,4) = sgn*3.0_dbl/4.0_dbl
   case('G'); SG%SYM_c(j-1,4) = sgn*5.0_dbl/6.0_dbl
   case('O'); SG%SYM_c(j-1,4) = 0.0_dbl
   case('X'); SG%SYM_c(j-1,4) = -sgn*3.0_dbl/8.0_dbl
   case('Y'); SG%SYM_c(j-1,4) = -sgn/4.0_dbl
   case('Z'); SG%SYM_c(j-1,4) = -sgn/8.0_dbl
  end select
 end do
end subroutine
!     
! ###################################################################
! 
!  subroutine MakeGenerators
!
!  Author: Marc De Graef
!  
!  Description: construct all generator matrices for the space group
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine MakeGenerators

use local
use crystalvars
use math

IMPLICIT NONE

integer(kind=irg),parameter :: QQ=48 
integer(kind=irg)           :: i,j,k,l,iset
real(kind=dbl)              :: SYM_d(4,4),SYM_e(4,4)
real(kind=dbl), parameter   :: eps=0.0005_dbl
character(1)                :: t(4)
character(40)               :: genst

! fill in the space group name 
 SG%SYM_name = SYM_SGname(cell%SYM_SGnum)
! initialize the encoded identity operator aOOO
 t(1)='a'
 do i=2,4
  t(i)='O'
 end do
! compute its matrix
 call SYM_fillgen(t,1)
! and put it first in the list of matrices
 do i=1,4
  do j=1,4
   SG%SYM_data(1,i,j) = SG%SYM_c(i,j)
  end do
 end do 
! get the space group generator string 
 genst = SYM_GL(cell%SYM_SGnum)
! initialize the number of generators 
 SG%SYM_GENnum = ichar(genst(2:2))-QQ
! create the generator matrices 
 do i=2,2+SG%SYM_GENnum - 1
   do  k=1,4
     l=2+4*(i-2)+k
     t(k) = genst(l:l)
   end do
   call SYM_fillgen(t,1)
   do k=1,4
    do l=1,4
     SG%SYM_data(i,k,l) = SG%SYM_c(k,l)
    end do
   end do
 end do
! this is where we are in the generator string
 i=2+4*SG%SYM_GENnum+1
! if there is inversion symmetry, add the inversion to the generators 
 if (genst(1:1).eq.'1') then 
  SG%SYM_centrosym=.TRUE.
  t(1)='h'
  do k=2,4
   t(k)='O'
  end do
  call SYM_fillgen(t,1)
  do j=1,4
   do k=1,4
    SG%SYM_data(SG%SYM_GENnum+2,j,k) = SG%SYM_c(j,k)
   end do
  end do
  SG%SYM_GENnum = SG%SYM_GENnum+2
 else   
  SG%SYM_GENnum = SG%SYM_GENnum+1
 end if
! now check for special origin conditions (choices 1 and 2) 
 if (genst(i:i).ne.'0') then 
  if (cell%SYM_SGset.eq.0) then
   call GetSetting(iset)
   cell%SYM_SGset=iset
  end if
  if (cell%SYM_SGset.eq.2) then 
! second setting: apply translation transformation to generators
   t(1)='a'
   do k=2,4
    l=i+k-1
    t(k) = genst(l:l)
   end do
   do l=2,SG%SYM_GENnum 
! translate to first setting origin
    call SYM_fillgen(t,-1)
    do j=1,4
     do k=1,4
      SYM_d(j,k)=SG%SYM_data(l,j,k)
     end do
    end do
! apply generator
    SYM_e = matmul(SYM_d,SG%SYM_c)
! translate back to second setting origin
    call SYM_fillgen(t,1)
    SYM_d = matmul(SG%SYM_c,SYM_e)
! reduce the translations to the fundamental unit cell
    do  k=1,3
     if (abs(SYM_d(k,4)).lt.eps) SYM_d(k,4)=0.0_dbl
     if (SYM_d(k,4).lt.0.0_dbl) then 
      do while (SYM_d(k,4).lt.0.0_dbl) 
       SYM_d(k,4)=SYM_d(k,4)+1.0_dbl
      end do
     end if
     if (SYM_d(k,4).ge.1.0_dbl) then 
      do while (SYM_d(k,4).ge.1.0_dbl) 
       SYM_d(k,4)=SYM_d(k,4)-1.0_dbl
      end do
     end if
     if (abs(SYM_d(k,4)-1.0_dbl).lt.eps) SYM_d(k,4)=0.0_dbl
    end do
! and store the result in the SYM_data array
    do j=1,4
     do k=1,4
      SG%SYM_data(l,j,k)=SYM_d(j,k)
     end do
    end do
   end do
  end if ! if (SYM_SGset.eq.2)
 end if ! if (genst(i:i).ne.'0')
end subroutine
!
! ###################################################################
! 
!  subroutine matrixmult
!
!  Author: Marc De Graef
!  
!  Description: multiplies two 4x4 symmetry matrices and brings
!               the translation component back to the fundamental
!               unit cell.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine matrixmult(k1,k2)
   
use local

IMPLICIT NONE

integer(kind=irg)        :: i,j,k,k1,k2
real(kind=dbl),parameter :: eps=0.0005_dbl

 do i=1,4
  do j=1,4
   SG%SYM_c(i,j) = 0.0_dbl
   do k=1,4
    SG%SYM_c(i,j)=SG%SYM_c(i,j)+SG%SYM_data(k1,i,k)*SG%SYM_data(k2,k,j)
   end do
  end do
 end do
! bring the translational part of the matrix back to
! the first unit cell and correct possible rounding errors
 do  k=1,3
  if (abs(SG%SYM_c(k,4)).lt.eps) SG%SYM_c(k,4)=0.0_dbl
  if (SG%SYM_c(k,4).lt.0.0_dbl) then 
   do while (SG%SYM_c(k,4).lt.0.0_dbl) 
    SG%SYM_c(k,4)=SG%SYM_c(k,4)+1.0_dbl
   end do
  end if
  if (SG%SYM_c(k,4).gt.1.0_dbl) then 
   do while (SG%SYM_c(k,4).gt.1.0_dbl) 
    SG%SYM_c(k,4)=SG%SYM_c(k,4)-1.0_dbl
   end do
  end if
  if (abs(SG%SYM_c(k,4)-1.0_dbl).lt.eps) SG%SYM_c(k,4)=0.0_dbl
 end do
end subroutine 
!
! ###################################################################
! 
!  function  isitnew
!
!  Author: Marc De Graef
!  
!  Description: make sure this is a new symmetry operator
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
logical function isitnew(nsym)

use local

IMPLICIT NONE

integer(kind=irg)        :: i,j,k,n,nsym
real(kind=dbl),parameter :: eps=0.0005_dbl

 k=0
 n=0
 do while ((k.le.nsym).and.(n.ne.12))
  n=0
  k=k+1
  do i=1,3
   do j=1,4
    if (abs(SG%SYM_c(i,j)- SG%SYM_data(k,i,j)).lt.eps) n=n+1
   end do
  end do
 end do
 if (n.ne.12) then 
  isitnew=.TRUE.
 else
  isitnew=.FALSE.
 end if
end function
!
! ###################################################################
! 
!  subroutine GenerateSymmetry
!
!  Author: Marc De Graef
!  
!  Description: compute all symmetry operators and store them 
!               in SG%SYM_data.
!               These routines are based on a program written 
!               by G. Ceder (MIT).
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine GenerateSymmetry(dopg)

use local
use crystalvars

IMPLICIT NONE

integer(kind=irg)    :: i,j,k,nsym,k1,k2,l1,l2
logical              :: dopg
real(kind=dbl)       :: q,sm

! create the space group generator matrices
 call MakeGenerators
 nsym = SG%SYM_GENnum
! generate new elements from the squares of the generators 
 do k=1,SG%SYM_GENnum 
  call matrixmult(k,k)
  if (isitnew(nsym).eqv..TRUE.) then 
   nsym=nsym+1
   do i=1,4
    do j=1,4
     SG%SYM_data(nsym,i,j) = SG%SYM_c(i,j)
    end do
   end do
  end if
 end do
! generate the remainder of the factorgroup
 k1=1
 do while (k1.le.nsym) 
  k2=k1+1
  do while (k2.le.nsym)
   call matrixmult(k2,k1)
   if (isitnew(nsym).eqv..TRUE.) then 
    nsym=nsym+1
    do i=1,4
     do j=1,4
      SG%SYM_data(nsym,i,j) = SG%SYM_c(i,j)
     end do
    end do
    if (nsym.ge.192) then 
     k2 = nsym
     k1 = nsym
    end if
   end if
   k2=k2+1
  end do
  k1=k1+1
 end do
 SG%SYM_MATnum = nsym
! reduce the translation operators to the fundamental unit cell
 do i=1,SG%SYM_MATnum
  do j=1,3
   SG%SYM_data(i,j,4)=mod( SG%SYM_data(i,j,4),1.0_dbl)
  end do
 end do
 if (dopg) then
! tag the point symmetry operators
! this is used to determine families of directions;
! for planes we must determine the transformed point symmetry
! operators SYM_recip() (this requires the metric tensors)
  SG%SYM_NUMpt=0
  do i=1,SG%SYM_MATnum 
   sm=SG%SYM_data(i,1,4)**2+SG%SYM_data(i,2,4)**2+SG%SYM_data(i,3,4)**2
   if (sm.lt.0.1_dbl) then
    SG%SYM_NUMpt=SG%SYM_NUMpt+1
! direct space point group symmetry elements
    do j=1,3
     do k=1,3 
      SG%SYM_direc(SG%SYM_NUMpt,j,k)=SG%SYM_data(i,j,k)
     end do
    end do
! reciprocal space point group symmetry elements
    do j=1,3
     do k=1,3 
      q=0.0_dbl
      do l1=1,3
       do l2=1,3
        q=q+cell%dmt(j,l1)*SG%SYM_data(i,l1,l2)*cell%rmt(l2,k)
       end do
      end do
      SG%SYM_recip(SG%SYM_NUMpt,j,k)=q
     end do
    end do
   end if  ! (sm.lt.0.1)
  end do
 end if ! if (dopg.eq..TRUE.)
! this completes generation of the factor group 
end subroutine
!
! ###################################################################
! 
!  subroutine CalcFamily
!
!  Author: Marc De Graef
!  
!  Description: compute the indices of equivalent planes/directions
!
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine CalcFamily(ind,num,space)
        
use local

IMPLICIT NONE

integer(kind=irg)          :: m,i,j,num,ind(3)
real(kind=sgl)             :: h,k,l,ih,ik,il,idiff
logical                    :: newpoint
character(1)               :: space
real,parameter             :: eps=0.0001

intent(IN)       :: ind,space
intent(OUT)      :: num

! first take the identity
 j=1
 itmp(j,1:3)=ind(1:3)
 h=float(ind(1))
 k=float(ind(2))
 l=float(ind(3))
! multiply with all point group elements
 do i=2,SG%SYM_NUMpt 
  if (space.eq.'d') then
   ih=SG%SYM_direc(i,1,1)*h+SG%SYM_direc(i,1,2)*k+SG%SYM_direc(i,1,3)*l
   ik=SG%SYM_direc(i,2,1)*h+SG%SYM_direc(i,2,2)*k+SG%SYM_direc(i,2,3)*l
   il=SG%SYM_direc(i,3,1)*h+SG%SYM_direc(i,3,2)*k+SG%SYM_direc(i,3,3)*l
  else
   ih=SG%SYM_recip(i,1,1)*h+SG%SYM_recip(i,1,2)*k+SG%SYM_recip(i,1,3)*l
   ik=SG%SYM_recip(i,2,1)*h+SG%SYM_recip(i,2,2)*k+SG%SYM_recip(i,2,3)*l
   il=SG%SYM_recip(i,3,1)*h+SG%SYM_recip(i,3,2)*k+SG%SYM_recip(i,3,3)*l
  end if
! is this a new point ?
  newpoint=.TRUE.
  do m=1,j+1
   idiff=(itmp(m,1)-ih)**2+(itmp(m,2)-ik)**2+(itmp(m,3)-il)**2
   if (idiff.lt.eps) newpoint=.FALSE.
  end do
  if (newpoint) then 
   j=j+1
   itmp(j,1)=nint(ih)
   itmp(j,2)=nint(ik)
   itmp(j,3)=nint(il)
  endif
 end do 
 num=j
end subroutine
!
! ###################################################################
! 
!  subroutine CalcOrbit
!
!  Author: Marc De Graef
!  
!  Description: computes the orbit of a given point (i.e., all
!               equivalent points in the fundamental unit cell)
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine CalcOrbit(m,n,ctmp)

use local
use crystalvars

IMPLICIT NONE

real(kind=dbl)            :: ctmp(192,3),r(3),s(3),diff
real(kind=dbl), parameter :: eps = 1.0D-4
integer(kind=irg)         :: n,i,j,k,mm,m
logical                   :: new

intent(IN)                :: m
intent(OUT)               :: n,ctmp

! get the atom coordinates
! and store them in the temporary array
 n=1
 do i=1,3
  r(i)=cell%ATOM_pos(m,i)
  ctmp(n,i)=r(i)
 end do
! get all the equivalent atom positions
 do i=2,SG%SYM_MATnum
  do j=1,3
   s(j)=SG%SYM_data(i,j,4)
   do k=1,3
    s(j)=s(j)+SG%SYM_data(i,j,k)*r(k)
   end do
  end do
! reduce to the fundamental unit cell if necessary
  if (SG%SYM_reduce) then
   do j=1,3
    s(j) = mod(s(j)+100.0_dbl,1.0_dbl)
   end do
  end if
! is this a new point ?
  new = .TRUE.
  do mm=1,n
   diff=0.0_dbl
   do j=1,3
    diff=diff+abs(ctmp(mm,j)-s(j))
   end do
   if (diff.lt.eps) then
     new = .FALSE.
   end if
  end do 
! yes, it is a new point
  if (new.eqv..TRUE.) then
   n=n+1
   do j=1,3
    ctmp(n,j)=s(j)
   end do
  end if
 end do
end subroutine
!
! ###################################################################
! 
!  subroutine CalcStar 
!
!  Author: Marc De Graef
!  
!  Description: computes the star of a given reciprocal vector
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   10/3/99 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine CalcStar(kk,n,stmp,space)

use local

IMPLICIT NONE

real(kind=dbl)           :: stmp(48,3),r(3),s(3),kk(3),diff
real(kind=dbl),parameter :: eps=1.0D-4
character(1)             :: space
logical                  :: new
integer(kind=irg)        :: n,i,j,k,mm

intent(IN)               :: kk,space
intent(OUT)              :: n,stmp

 n=1
 r=kk
 stmp(n,1:3)=r(1:3)
! get all the equivalent reciprocal/direct space vectors
 do i=2,SG%SYM_NUMpt 
  do j=1,3
   s(j)=0.0_dbl
   do k=1,3
    if (space.eq.'r') then 
     s(j)=s(j)+SG%SYM_recip(i,j,k)*r(k)
    else
     s(j)=s(j)+SG%SYM_direc(i,j,k)*r(k)
    end if
   end do
  end do
! is this a new point ?
  new = .TRUE.
  do mm=1,n
   diff=0.0_dbl
   do j=1,3
    diff=diff+abs(stmp(mm,j)-s(j))
   end do
   if (diff.le.eps) then
     new  = .FALSE.
   endif
  end do
! yes, it is a new point
  if (new.eqv..TRUE.) then
   n=n+1
   stmp(n,1:3)=s(1:3)
  end if
 end do
end subroutine
!
! ###################################################################
! 
!  subroutine CalcPositions
!
!  Author: Marc De Graef
!  
!  Description: Compute all atom positions in the fundamental unit
!               cell and translate to neighbouring cells if needed
!               (used for structure drawings and structure factor
!               computations)
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine CalcPositions(switch)

use local
use crystalvars
use io
use error
use crystal

IMPLICIT NONE

logical           :: inside
integer(kind=irg) :: i,j,k,l,mm,icnt,celln(3),ncells,n,kk,ier
character(*),parameter :: var1 = 'CalcPositions: unable to allocate memory for array '
real(kind=dbl)    :: ctmp(192,3),ff(3),sh(3)
real(kind=sgl)    :: r(3),g(3)
character(1)      :: switch

intent(IN)        :: switch

! make sure all coordinates are reduced to the fundamental unit cell
 SG%SYM_reduce=.TRUE.
! multiple cells ?
 if (switch.eq.'m') then 
  mess = 'Number of unit cells in a, b and c direction ?: '; call GetInt(3)
  do j=1,3
   celln(j) = io_int(j)
   sh(j) = 0.5_dbl*celln(j)+1.0_dbl
  end do
  ncells = celln(1)*celln(2)*celln(3)
 else
! no, just one cell
  do j=1,3
   celln(j)=0
   sh(j)=1.0_dbl
  end do
  ncells = 1
 end if
! main loop
! first allocate the apos variable (contains CARTESIAN coordinates
! if switch is 'm', crystal coordinates otherwise)
 if (allocated(apos)) deallocate(apos)
 allocate (apos(cell%ATOM_ntype, ncells * SG%SYM_MATnum, 3),stat=ier)
 if (ier.ne.0) call FatalError(var1,'apos')
!
 do i=1,cell%ATOM_ntype
! for each atom in the asymmetric unit
  call CalcOrbit(i,n,ctmp)
  numat(i)=n
  icnt=1
! replicate in all cells
  do j=1,celln(1)+1
   ff(1)=float(j)
   do k=1,celln(2)+1
    ff(2)=float(k)
    do l=1,celln(3)+1
     ff(3)=float(l)
     do kk=1,numat(i)
      do mm=1,3
       r(mm)=ctmp(kk,mm)+ff(mm)-sh(mm)
      end do 
      if (switch.eq.'m') then
! make sure the atom is actually inside the block of unit
! cells, or on one of the edges/faces/corners
       inside=.TRUE.
       do mm=1,3
        if ((r(mm)+sh(mm)).gt.(celln(mm)+1.0)) inside=.FALSE.
       end do
       if (inside) then
        call TransSpace(r,g,'d','c')
        do mm=1,3
         apos(i,icnt,mm)=g(mm)
        end do
        icnt=icnt+1
       end if
      else
! prepare for structure factor computation
       do mm=1,3
        apos(i,icnt,mm)=r(mm)
       end do
       icnt=icnt+1
      end if
     end do 
    end do 
   end do 
  end do 
  numat(i)=icnt-1
 end do
end subroutine
! ###################################################################
!  subroutine GetSetting
!
!                                    created: 10/13/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: for the space groups with a second origin setting
!               this routine asks which of the two settings to use
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
subroutine GetSetting(iset)

use local
use io
use crystalvars

IMPLICIT NONE

integer(kind=irg)      :: i,iset,isg
! There are 24 space groups with two origin choices.
! The symmetry of both sites is stored in the array
! sitesym; the space group numbers are stored
! in tworig
integer(kind=irg),parameter  :: tworig(24)=(/48,50,59,68,70,85,86,88,125,126,129,130,133,134,137,138,&
                                            141,142,201,203,222,224,227,228/)
character(7),parameter :: sitesym(48) = (/ '222    ',' -1    ','222/n  ',' -1    ','mm2/n  ',' -1    ', &
                                           '222    ',' -1    ','222    ',' -1    ','-4     ',' -1    ', &
                                           '-4     ',' -1    ','-4     ',' -1    ','422    ','2/m    ', &
                                           '422/n  ',' -1    ','-4m2   ','2/m    ','-4/ncn ',' -1    ', &
                                           '-4121/c',' -1    ','-42m   ','2/m    ','-4m2/n ',' -1    ', &
                                           '-4cg   ','2/m    ','-4m2   ','2/m    ','-4c21  ',' -1    ', &
                                           '23     ',' -3    ','23     ',' -3    ','432    ',' -3    ', &
                                           '-43m   ','-3m    ','-43m   ','-3m    ','23     ',' -3    '/)

intent(OUT)            :: iset

 isg = 0
 do i=1,24
  if (tworig(i).eq.cell%SYM_SGnum) isg=i
 end do
 if (isg.ne.0) then 
  mess = '---------------------------------------------'; call Message("(A)")
  mess = 'This space group has two origin settings.'; call Message("(A)")
  mess = 'The first setting has site symmetry    : '//sitesym(2*isg-1); call Message("(A)")
  mess = 'the second setting has site symmetry   : '//sitesym(2*isg); call Message("(A)")
  mess = 'Which setting do you wish to use (1/2) : '; call GetInt(1)
  iset = io_int(1)
  mess = '---------------------------------------------'; call MEssage("(A)")
 end if
end subroutine
!
! ###################################################################
! 
!  subroutine GetSpaceGroup
!
!                                    created: 10/13/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: asks for the space group (lists appropriate groups)
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
!
subroutine GetSpaceGroup

use local
use io
use crystalvars

IMPLICIT NONE

integer(kind=irg)  :: sgmin,sgmax,i,j,TRIG(7)
logical            :: skip

!
! the following space groups have a hexagonal and rhombohedral
! setting;  the generator matrices are different, and so are
! their entries in the SYM_GL array.
!  hexagonal setting 146: R 3           231
!  hexagonal setting 148: R -3          232
!  hexagonal setting 155: R 3 2         233
!  hexagonal setting 160: R 3 m         234
!  hexagonal setting 161: R 3 c         235
!  hexagonal setting 166: R -3 m        236
!  hexagonal setting 167: R -3 c        237
!
 TRIG = (/ 146,148,155,160,161,166,167 /)
 skip = .FALSE.
 select case (cell%xtal_system)
 case (1); sgmin = 195; sgmax = 230
 case (2); sgmin =  75; sgmax = 142
 case (3); sgmin =  16; sgmax =  74
 case (4); sgmin = 168; sgmax = 194
 case (5); if (cell%SYM_second) then
             mess = 'The space groups below correspond to the '; call Message("(A)")
             mess = 'second (rhombohedral) setting.'; call Message("(A/)")
             mess = 'Please select one of the space groups.'; call Message("(A/)")
             do i=1,7
              if ((mod(i,4).eq.0).or.(i.eq.7)) then
                write (*,"(1x,i3,':',A11,5x)") TRIG(i),SYM_SGname(TRIG(i))
              else
                write (*,"(1x,i3,':',A11,5x,$)") TRIG(i),SYM_SGname(TRIG(i))
              end if
             end do 
             mess = ' -------------------------- '; call Message("(A)")
             mess = ' Enter space group number : '; call GetInt(1)
             cell%SYM_SGnum = io_int(1)

! check for rhombohedral settings of rhombohedral space groups
             if (cell%SYM_second) then
               if (cell%SYM_SGnum.eq.146) cell%SYM_SGnum=231
               if (cell%SYM_SGnum.eq.148) cell%SYM_SGnum=232
               if (cell%SYM_SGnum.eq.155) cell%SYM_SGnum=233
               if (cell%SYM_SGnum.eq.160) cell%SYM_SGnum=234
               if (cell%SYM_SGnum.eq.161) cell%SYM_SGnum=235
               if (cell%SYM_SGnum.eq.166) cell%SYM_SGnum=236
               if (cell%SYM_SGnum.eq.167) cell%SYM_SGnum=237
             endif
             skip = .TRUE.
           else 
            sgmin = 143
            sgmax = 167
           end if
 case (6); sgmin =   3; sgmax =  15
 case (7); sgmin =   1; sgmax =   2
 end select
! print out all the relevant space group names and numbers        
 if (skip.eqv..FALSE.) then
  do i=sgmin,sgmax
   j=i-sgmin+1
   if ((mod(j,4).eq.0).or.(i.eq.sgmax)) then
    write (*,"(1x,i3,':',A11,5x)") i,SYM_SGname(i)
   else
    write (*,"(1x,i3,':',A11,5x,$)") i,SYM_SGname(i)
   end if
  end do
  cell%SYM_SGnum = sgmin-1
  do while ((cell%SYM_SGnum.lt.sgmin).or.(cell%SYM_SGnum.gt.sgmax)) 
   mess = ' -------------------------- '; call Message("(A)")
   mess = ' Enter space group number : '; call GetInt(1)
   cell%SYM_SGnum = io_int(1)
   if ((cell%SYM_SGnum.lt.sgmin).or.(cell%SYM_SGnum.gt.sgmax)) then
    mess = 'Error in space group number '; call Message("(A)")
    mess = 'Crystal system / space group mismatch '; call Message("(A)")
   end if
  end do
 end if 
end subroutine
! ###################################################################
! 
!  subroutine GetOrder  
! 
!                                    created: 5/24/01 
!  Author: Marc De Graef
!  
!  Description: determine the order of the subfamily of a reciprocal 
!               lattice family belonging to a zone
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   5/24/01 MDG 1.0 original
! 11/27/01  MDG 1.1 added kind support
! ###################################################################
subroutine GetOrder(k,il,num,jcnt)

use local
use symmetryvars   ! uses itmp variable

IMPLICIT NONE

integer(kind=irg)        :: i,jcnt,il(48),num
real(kind=sgl)           :: k(3),gn(3)
real(kind=sgl),parameter :: eps=1.0E-5

intent(IN)     :: k,num
intent(OUT)    :: il,jcnt

 jcnt = 0
 do i=1,num
  gn(1:3) = float(itmp(i,1:3))
  if (abs(sum(gn*k)).lt.eps) then
    jcnt = jcnt+1
    il(jcnt) = i
  end if
 end do
end subroutine 
! ###################################################################
! 
!  subroutine ShortestG
!
!                                    created: 10/13/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: determine the pair of shortest reciprocal lattice
!               vectors for a given zone axis
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/24/01 MDG 2.0 f90 version
!   9/25/01 MDG 2.1 added 2D point group stuff
! 11/27/01  MDG 2.2 added kind support
! ###################################################################
!
subroutine ShortestG(k,gone,gtwo,isym)

use local
use error
use crystal
use crystalvars
use constants

IMPLICIT NONE

integer(kind=irg)      :: ga(3),gb(3),nzero(3),k(3),gone(3),gtwo(3),u,v,w,snz,ml(4),seo,igsave(3)
integer(kind=irg)      :: ima,imb,ina,inb,el(6),denom,minsum,inm,i,isym,il(48),jcnt,num
real(kind=sgl)         :: fel(6),fit,gsave(3)
integer(kind=irg),allocatable :: ifit(:,:,:,:)

INTENT(IN)    :: k
INTENT(OUT)   :: gone, gtwo
INTENT(INOUT) :: isym

 u = k(1)
 v = k(2)
 w = k(3)
! determine two arbitrary vectors normal to k 
! first count the zeroes in k
 nzero = (/0,0,0/)
 where (k.eq.0) nzero = 1
 snz = sum(nzero)
 if (snz.eq.0) then  ! make sure ga x gb is parallel to k
   ga = (/v,-u,0/)
   gb = (/w,0,-u/)
 else
  select case (snz)
  case(1);  ga = nzero; gb = 1 - nzero
            if (nzero(1).eq.1) gb = gb * (/0,w,-v/)
            if (nzero(2).eq.1) gb = gb * (/w,0,-u/)
            if (nzero(3).eq.1) gb = gb * (/v,-u,0/)
  case(2);  if ((nzero(1).eq.1).and.(nzero(2).eq.1)) then
             ga = (/1,0,0/); gb = (/0,1,0/)
            endif
            if ((nzero(1).eq.1).and.(nzero(3).eq.1)) then
             ga = (/0,0,1/); gb = (/1,0,0/)
            endif
            if ((nzero(2).eq.1).and.(nzero(3).eq.1)) then
             ga = (/0,1,0/); gb = (/0,0,1/)
            endif
  case(3); call ErrorMess(5); stop
  end select
 end if 
! next look for 4 integer numbers which transform ga and gb
! simultaneously into two shortest possible vectors of the same zone;
! a range of -5:5 should be sufficient as a search space
! 
! If g1 and g2 are those shortest vectors, then we have
! 
!    ga  =  na*g1 + ma*g2
!    gb  =  nb*g1 + mb*g2
! 
! Inversion of this relation gives
! 
!    g1  =  (mb*ga - ma*gb)/D
!    g2  =  (-nb*ga + na*gb)/D
! 
! with D = na*mb - nb*ma
! 
! The procedure below searches for the combination of 
! (ma,mb,na,nb) which simultaneously minimizes the 
! length of g1 and g2, and makes sure that both g1 and g2
! are integer linear combinations of the reciprocal basis
! vectors.
! 
 inm = 5
 allocate(ifit(-inm:inm,-inm:inm,-inm:inm,-inm:inm))
 do ima=-inm,inm
  do imb=-inm,inm
   do ina=-inm,inm
    do inb=-inm,inm
     el(1) = imb*ga(1)-ima*gb(1) 
     el(2) = imb*ga(2)-ima*gb(2) 
     el(3) = imb*ga(3)-ima*gb(3) 
     el(4) = ina*gb(1)-inb*ga(1) 
     el(5) = ina*gb(2)-inb*ga(2) 
     el(6) = ina*gb(3)-inb*ga(3) 
     denom = ina*imb-inb*ima
     ifit(ima,imb,ina,inb)=100
     if (denom.ne.0) then
      fel = float(el)/float(denom)
      fit = sum(abs(float(int(fel))-fel))
! here is where we only keep the integer combinations
      if (fit.eq.0.0) then
        gone(1:3) = int(fel(1:3))
        gtwo(1:3) = int(fel(4:6))
! keep the sum of the squares of the lengths 
       ifit(ima,imb,ina,inb)=sum(gone**2)+sum(gtwo**2) 
      end if
     end if
    end do
   end do
  end do
 end do
 minsum = 50
! look for the minimum of ifit with the smallest and most
! positive coefficients; store them in ml
! [minloc does not work here because there may be multiple minima]
 do ima=-inm,inm
  do imb=-inm,inm
   do ina=-inm,inm
    do inb=-inm,inm
     if (ifit(ima,imb,ina,inb).le.minsum) then
      minsum = ifit(ima,imb,ina,inb)
      ml(1) = ima
      ml(2) = imb
      ml(3) = ina
      ml(4) = inb
     end if
    end do
   end do
  end do
 end do
 deallocate(ifit)
! finally transform ga and gb into g1 and g2 
 gone = (ml(2)*ga-ml(1)*gb)/(ml(3)*ml(2)-ml(4)*ml(1))
 gtwo = (ml(3)*gb-ml(4)*ga)/(ml(3)*ml(2)-ml(4)*ml(1))
! next rank these two vectors so that their cross product is along +k
 call CalcCross(float(gone),float(gtwo),gsave,'r','r',0)
 fit = CalcDot(gsave,float(k),'r')
 if (fit.lt.0.0) then
  igsave = gone
  gone = gtwo
  gtwo = igsave
 end if

! finally, if isym.ne.0 make sure that the selection of the 
! basis vectors for the 3-fold and 6-fold 2D point groups is
! correctly done.
!
! For isym < 7:   90 degrees between vectors (should not be a problem)
! For isym = 7:  120 degrees between vectors (should be ok too)
! distinguish between 3 coming from a cubic group
! vs. the same symmetries originating from a hexagonal setting
 if ((isym.eq.7).and.(cell%gamma.eq.120.0)) then
   isym=11
   fit = CalcAngle(float(gone),float(gtwo),'r')*180.0/cPi 
   if (abs(fit-120.0).lt.1.0) then
     gtwo=gone+gtwo
   end if
 end if
! For isym = 8:  here we should distinguish between the settings 3m1 and 31m !!!
!                The angle should always be 120 degrees, so we must check that
!                this is the case for the selected gone and gtwo.
 if ((isym.eq.8).and.(cell%gamma.eq.120.0)) then
   isym=12
   fit = CalcAngle(float(gone),float(gtwo),'r')*180.0/cPi 
   if (abs(fit-120.0).lt.1.0) then
     gtwo=gone+gtwo
   end if
 end if
!
 if (isym.eq.8) then
   fit = CalcAngle(float(gone),float(gtwo),'r')*180.0/cPi 
   if (abs(fit-120.0).gt.1.0) then
     gtwo=gtwo-gone
   end if
! we assume it is the 31m setting;  if the order of gone is 6, then that is not true
   call CalcFamily(gone,num,'r')
   call GetOrder(float(k),il,num,jcnt)
   if (jcnt.eq.6) then  ! it is the 3m1 setting
     isym = 13
   end if
 end if
! it could the 3m1 setting for the 3m hexagonal case
 if (isym.eq.12) then
! we assume it is the 31m setting;  if the order of gone is 6, then that is not true
   call CalcFamily(gone,num,'r')
   call GetOrder(float(k),il,num,jcnt)
   if (jcnt.eq.6) then  ! it is the 3m1 setting
     isym = 14
   end if
 end if
! For isym = 9 or 10:   60 degrees between vectors (may not be the case)
 if ((isym.eq.9).or.(isym.eq.10)) then
   fit = CalcAngle(float(gone),float(gtwo),'r')*180.0/cPi 
   if (abs(fit-120.0).lt.1.0) then
     gtwo=gone+gtwo
   end if
 end if

end subroutine
!
! ###################################################################
! 
!  function IsGAllowed
!
!                                    created:  5/24/01
!  Author: Marc De Graef
!  
!  Description: logical function to determine whether or not a 
!               reflection is allowed by lattice centering
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   5/24/01 MDG 1.0 f90 version
! 11/27/01  MDG 2.1 added kind support
! ###################################################################
logical function IsGAllowed(g)

use local
use symmetryvars

IMPLICIT NONE

integer(kind=irg)       :: g(3), seo
character(1)            :: lc

intent(IN)              :: g

! Determine whether or not this vector is
! actually allowed by the lattice centering
 lc(1:1) =  SG%SYM_name(2:2)
 IsGAllowed = .TRUE.
 select case (lc)
  case ('P'); ! all reflections allowed for a primitive lattice
  case ('F'); seo = sum(mod(g+100,2)); if ((seo.eq.1).or.(seo.eq.2)) IsGAllowed = .FALSE.
  case ('I'); seo = mod(sum(g)+100,2); if (seo.eq.1) IsGAllowed = .FALSE.
  case ('A'); seo = mod(g(2)+g(3)+100,2); if (seo.eq.1) IsGAllowed = .FALSE.
  case ('B'); seo = mod(g(1)+g(3)+100,2); if (seo.eq.1) IsGAllowed = .FALSE.
  case ('C'); seo = mod(g(1)+g(2)+100,2); if (seo.eq.1) IsGAllowed = .FALSE.
  case ('R'); if (hexset) then
               seo = mod(-g(1)+g(2)+g(3)+90,3); if (seo.ne.0) IsGAllowed = .FALSE.
              endif ! otherwise reflections are all allowed
 end select
end function
!
! ###################################################################
! 
!  subroutine BFsymmetry
!
!                                    created:  5/24/01
!  Author: Marc De Graef
!  
!  Description: routine to determine the symmetry of a bright field
!               bend center [uses Table 4 in BESR paper]
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   6/22/01 MDG 1.0 f90 version
! 11/27/01  MDG 1.1 added kind support
! ###################################################################
subroutine BFsymmetry(uvw,j,isym,ir)

use local

IMPLICIT NONE

integer(kind=irg)     :: orderPG, Lauenum, ir, ng, uvw(3), isym,j
real(kind=dbl)        :: kstar(48,3)

intent(IN)  :: uvw

 orderPG = SG%SYM_numpt
 Lauenum = PGLaueinv(j)
 call CalcStar(dble(uvw),ng,kstar,'d')
 ir = orderPG/ng
! intercept the special cases where more than one symmetry of the same order can
! occur in a given Laue group;  e.g.  Laue group 2/m has two bright field symmetries
! of order 2:  2 and m
! 
! This happens for the following Laue groups:
!   2/m     [010] -> 2    [u0w] -> m
!   -3m     [11.0] -> 2   [u-u.w] -> m
! 
! The normal conversion from the reduced order ir to the actual Bright Field
! symmetry uses the PGTWDinverse array to determine the 2D symmetry
! 
 isym = PGTWDinverse(ir,Lauenum)
 if ((Lauenum.eq.2).and.(ir.eq.2)) then   ! this deals with Laue Group 2/m
  if (uvw(2).eq.0) isym=3
 end if
 if ((Lauenum.eq.7).and.(ir.eq.2)) then   ! and this covers -3m
  if (uvw(1).eq.-uvw(2)) isym=3
 end if
end subroutine

end module

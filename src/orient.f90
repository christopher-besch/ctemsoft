!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:  orient.f90                                                         !
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
!  FILE: "orient.f90"
!                                    created: 10/13/98 {9:29:46 AM} 
!                               
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: stereographic projection of a crystal orientation relation
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! contains errors that need to be corrected
! ###################################################################
program orient

use local
use crystalvars
use crystal
use symmetryvars
use graphics
use files
use postscript
use io
use math

character(1)      :: sp
logical           :: nn,topbot
type(unitcell)    :: cellA,cellB
type(orientation) :: orel
real(kind=sgl)    :: rr(3),gg(3),g(3),r(3),M(3,3),negthresh,p(3),Ep(3,3),E(3,3),TT(3,3),q(3)
real(kind=dbl)    :: dE(3,3),dgg(3)
integer(kind=irg) :: h,k,l,cr,hkl(3),iview(3)

 progname = 'orient.f90'
 progdesc = 'Stereographic projection of orientation relation'
 call CTEMsoft

 inm=2
 SG % SYM_reduce=.TRUE.
 topbot=.FALSE.
! 20cm radius projection circle, centered on page [inches]
 CRad = 3.937
 CX = 3.25
 CY = 3.5
 negthresh=-0.0001
! read crystal A information
 call CrystalData
! store Crystal A matrices
 cella = cell
! read crystal B information
 call CrystalData
! store Crystal B matrices 
 cellb = cell
! get orientation relation
 call GetOR(orel)
! compute E matrix  [page 74]
 cell = cella
 call TransSpace(orel % gA,r,'r','d')
 call NormVec(r,'d')
 call NormVec(orel % tA,'d')
 call CalcCross(orel % tA,r,p,'d','d',0)
 call NormVec(p,'d')
 do i=1,3
  E(1,i)=r(i)
  E(2,i)=p(i)
  E(3,i)=orel % tA(i)
 end do
 call mInvert(dble(E),dE,.FALSE.)
 E = dE
 mess = 'Transformation matrix E'; call Message("(A)")
 do i=1,3
  write (*,*) (E(i,j),j=1,3)
 end do
! compute E-prime matrix 
 cell = cellb
 call TransSpace(orel % gB,r,'r','d')
 call NormVec(r,'d')
 call NormVec(orel % tB,'d')
 call CalcCross(orel % tB,r,p,'d','d',0)
 call NormVec(p,'d')
 do i=1,3
  Ep(1,i)=r(i)
  Ep(2,i)=p(i)
  Ep(3,i)=orel % tB(i)
 end do
 mess ='Transformation matrix E-prime'; call Message("(A)")
 do i=1,3
  write (*,*) (Ep(i,j),j=1,3)
 end do
! and multiply both matrices to get transformation matrix M
 TT = matmul(E,Ep)
 mess = 'Transformation matrix for orientation relation'; call Message("(A)")
 do i=1,3
  write (*,*) (TT(i,j),j=1,3)
 end do
 write (*,*) ' -- '
! real space or reciprocal space
 call GetDrawingSpace(sp)
! and from here one it is the same as a regular stereographic projection
! except that there are two sets of points to be drawn.
! viewing direction
 call GetViewingDirection(iview)
! create transformation matrix
 call ProjectionMatrix(iview,M)
! open PostScript file
 call PS_openfile
! write text and draw projection circle
 call DrawFrame(CX,CY,CRad,iview,sp,cella,cellb,orel)
! loop over all planes or directions
 do h=-inm,inm
  do k=-inm,inm
   do l=-inm,inm
    ih=h
    ik=k
    il=l
! skip the origin
    if ((ih**2+ik**2+il**2).ne.0) then
! reduce to smallest integers to avoid overlap
! of indices, such as (111) and (222)
     hkl(1)=ih
     hkl(2)=ik
     hkl(3)=il
     call IndexReduce(hkl)
! transform to cartesian coordinates
     g(1)=float(hkl(1))
     g(2)=float(hkl(2))
     g(3)=float(hkl(3))
     ih = hkl(1)
     ik = hkl(2)
     il = hkl(3)
! crystal A
     cell = cella
     call TransSpace(g,r,sp,'c')
     call NormVec(r,'c')
! apply viewing tansformation
     rr = matmul(M,r)
! compute stereographic projection coordinates
     xst=CX+CRad*rr(1)/(1.0+abs(rr(3)))
     yst=CY+CRad*rr(2)/(1.0+abs(rr(3)))
     cr=1
     if (rr(3).gt.negthresh) then
      call PS_filledcircle(xst,yst,0.015/PS % psscale,0.0)
      nn = .TRUE.
      call DumpIndices(sp,ih,ik,il,cr,xst,yst,nn)
     else if (topbot) then
      call PS_circle(xst,yst,0.035/PS % psscale)
      nn = .FALSE.
      call DumpIndices(sp,ih,ik,il,cr,xst,yst,nn)
     end if
! crystal B
     cell = cellb
     call TransCoor(dble(g),dgg,dble(TT),sp,'on')
     gg = sngl(dgg)
     call TransSpace(gg,r,sp,'c')
     call NormVec(r,'c')
! apply viewing tansformation
     rr = matmul(M,r)
! compute stereographic projection coordinates
     xst=CX+CRad*rr(1)/(1.0+abs(rr(3)))
     yst=CY+CRad*rr(2)/(1.0+abs(rr(3)))
     cr=2
     if (rr(3).gt.negthresh) then
      call PS_filledsquare(xst,yst,0.035/PS % psscale,0.0)
      nn = .TRUE.
      call DumpIndices(sp,ih,ik,il,cr,xst,yst,nn)
     else if (topbot) then
      call PS_square(xst,yst,0.050/PS % psscale)
      nn = .FALSE.
      call DumpIndices(sp,ih,ik,il,cr,xst,yst,nn)
     end if
    end if
   end do 
  end do 
 end do 
! close Postscript file
 call PS_closefile
end program
! ###################################################################
!
!  subroutine DrawFrame
!
!  Author: Marc De Graef
!
!  Description: draw a stereographic projection layout
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
subroutine DrawFrame(CX,CY,CRad,iview,sp,cella,cellb,orel)

use local
use io
use postscript
use crystalvars

type(unitcell)     :: cella,cellb
type(orientation)  :: orel
character(17)      :: str
character(12)      :: instr
integer            :: hkl(3),iview(3)
real               :: CX, CY, CRad
character(1)       :: sp

 call PS_newpage(.FALSE.,'Stereographic Projection')
 call PS_setlinewidth(0.016)
 call PS_circle(CX,CY,CRad)
 call PS_setlinewidth(0.004)
 call PS_line(CX-CRad,CY,CX+CRad,CY)
 call PS_line(CX,CY-CRad,CX,CY+CRad)
 call PS_setfont(PSfonts(2),0.08)
 call PS_text(CX-CRad-0.07,CY-0.025,'A')
 call PS_text(CX+CRad+0.03,CY-0.025,'B')
 call PS_text(CX-0.03,CY-CRad-0.09,'M''')
 call PS_text(CX-0.03,CY+CRad+0.07,'M"')
 call PS_setfont(PSfonts(2),0.12/PS % psscale)
 call PS_text(0.35,8.30,'Crystal A : '//cella % fname)
 call PS_filledcircle(0.0,8.30,0.015,0.0)
 call DumpIndices(sp,0,0,0,1,0.0,8.30,.TRUE.)
 call PS_setfont(PSfonts(2),0.12)
 call PS_text(0.35,8.10,'Crystal B : '//cellb % fname)
 call PS_filledsquare(0.0,8.10,0.035,0.0)
 call DumpIndices(sp,0,0,0,2,0.0,8.10,.TRUE.)
 call PS_setfont(PSfonts(2),0.12)
 call IndexString(instr,iview,'d')
 call PS_text(0.0,7.90,'Viewing Direction '//instr//' [A]')
 if (sp.eq.'d') then 
  str='direct space'
 else
  str='reciprocal space'
 endif
 call PS_text(0.0,7.70,'Projection of '//str)
!
 call PS_text(CX,8.20,'Orientation Relation ')
 do i=1,3 
  hkl(i)=int(orel % gA(i))
 end do
 call IndexString(instr,hkl,'r')
 call PS_text(CX,8.00,'\(hkl\) : ')
 call PS_text(CX+0.4,8.00,'A-'//instr)
 do i=1,3 
  hkl(i)=int(orel % gB(i))
 end do
 call IndexString(instr,hkl,'r')
 call PS_text(CX+0.9,8.00,'|| B-'//instr)
! Space=.True.
 do i=1,3 
  hkl(i)=int(orel % tA(i))
 end do 
 call IndexString(instr,hkl,'d')
 call PS_text(CX,7.80,'[uvw] : ')
 call PS_text(CX+0.4,7.80,'A-'//instr)
 do i=1,3 
  hkl(i)=int(orel % tB(i))
 end do
 call IndexString(instr,hkl,'d')
 call PS_text(CX+0.9,7.80,'|| B-'//instr)
end subroutine


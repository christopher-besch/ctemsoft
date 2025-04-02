!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft: drawcell.f90                                                                    !
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
!  FILE: "drawcell.f90"
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/5/99 {10:37:32 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: draw a unit cell
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
program drawcell 

use local
use io
use crystalvars
use crystal
use symmetryvars
use symmetry
use graphics
use postscript
use math
use constants
use files
        
integer(kind=irg),parameter    :: n=1000
character(1)                   :: sgn,ans,sp
character(3)                   :: acol(n)
logical                        :: first,again,nn,more,topbot,inside
real(kind=sgl)                 :: p(4),q(4),xmax,x(n),y(n),z(n),x1,y1,z1,asize(n),rr(4),gg(4),g(3), &
                                  r(3),cc,M(4,4),negthresh,VD,sc,diam
integer(kind=irg)              :: h,k,l,mi,kk,idx(n),iview(3),iform

 progname = 'drawcell.f90'
 progdesc='Draw one or more unit cells in perspective mode'
 call CTEMsoft
 
 SG % SYM_reduce=.TRUE.
 CX=7.0
 CY=7.0
! read crystal information
 call CrystalData
! real space or reciprocal space
 call GetDrawingSpace(sp)
! create all atoms
 call CalcPositions('m')
! Viewing Distance = 3 times CX = 6 times xmax
 mess = 'Viewing distance :'; call GetReal(1); VD = io_real(1)
! viewing direction
 call GetViewingDirection(iview)
! viewing direction
!       call ProjectionMode(iview)
! create transformation matrix
 call ComputeViewTrans(iview,M,VD)
! open PostScript file
 call PS_openfile
! write text and draw box
 call DrawFrame(iview,sp,CX,CY)
! draw unit cell outline first
! which radii should be used for the drawing ?
 mess = 'use ionic radii (1) or metallic (2) :'; call GetInt(1); iform = io_int(1)
! then get all atom coordinates
 icnt=0
 xmin=100.0
 xmax=-100.0
 ymin=100.0
 ymax=-100.0
 do i=1,cell % ATOM_ntype
  do j=1,numat(i)
  write (*,fmt="(3(f10.5,2x))") apos(i,j,1),apos(i,j,2),apos(i,j,3)
   p=(/sngl(apos(i,j,1)),sngl(apos(i,j,2)),sngl(apos(i,j,3)),1.0/)
   q = matmul(p,M)
   x1=VD*q(1)/q(4) 
   y1=VD*q(2)/q(4) 
   z1=VD*q(3)/q(4) 
   x1 = VD*x1/z1
   y1 = VD*y1/z1
   icnt=icnt+1
   x(icnt)=0.5*CX+2.5*x1
   y(icnt)=0.5*CY+2.5*y1
   if (x(icnt).lt.xmin) xmin=x(icnt)
   if (x(icnt).gt.xmax) xmax=x(icnt)
   if (y(icnt).lt.ymin) ymin=y(icnt)
   if (y(icnt).gt.ymax) ymax=y(icnt)
   z(icnt)=z1
   if (iform.eq.1) then 
    asize(icnt)=ATOM_SPradii(cell % ATOM_type(i))
   else
    asize(icnt)=ATOM_MTradii(cell % ATOM_type(i))
   endif
   acol(icnt)=ATOM_color(cell % ATOM_type(i))
  end do
 end do
! shift the drawing back to the center of the page
 xmid = (xmax+xmin)*0.5
 ymid = (ymax+ymin)*0.5
 shx = 0.5*CX - xmid
 shy = 0.5*CY - ymid
 do i=1,icnt
  x(i) = x(i) + shx
  y(i) = y(i) + shy
  write (*,fmt="(2(f10.5,2x))") x(i),y(i)
 end do
! rank according to distance from observer
! draw atoms/relpoints in reverse order (farthest first)
 call SPSORT(z,icnt,idx,-1,ier)
write(*,fmt="(I5)") icnt
 do i=1,icnt
  j=idx(i) 
  diam = asize(j)
  call PS_sphere(x(j),y(j),diam,acol(j))
  write (*,fmt="(3(f10.5,2x),A3)") x(j),y(j),diam,acol(j)
 end do 
! close PostScript file
 call PS_closefile
! formats
end program
! ###################################################################
!
!  subroutine DrawFrame
!
!  Author: Marc De Graef
!
!  Description: format the page 
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
subroutine DrawFrame(iview,sp,CX,CY)

use local
use io
use postscript

logical           :: SaveSpace
character(12)     :: instr
character(17)     :: str
character(1)      :: sp
integer(kind=irg) :: iview(3)

 if (sp.eq.'d') then 
  call PS_newpage(.FALSE.,'Crystal Structure Drawing')
 else
  call PS_newpage(.FALSE.,'Reciprocal Lattice Drawing')
 endif
 call PS_setlinewidth(0.012)
 call PS_drawrect(0.0,0.0,CX,CY)
 call PS_setlinewidth(0.008)
 call PS_setfont(PSfonts(2),0.12/PS % psscale)
 call PS_cellinfo(0.0,8.3)
 call IndexString(instr,iview,'d')
 call PS_text(CX*0.5,8.14,'Viewing Direction '//instr)
 if (sp.eq.'d') then 
  str='direct space'
 else
  str='reciprocal space'
 endif
 call PS_text(CX*0.5,8.00,'Drawing of '//str)
end subroutine

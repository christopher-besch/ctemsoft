!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:postscript.f90                                                       !
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
!  FILE: "postscript.f90"
!                                    created: 12/31/92 {12:47:01 PM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: 
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  12/31/92 MDG 1.0 original
!  1/8/98   MDG 2.0 added new routines
!  7/19/99  MDG 2.1 added sp command tp Postscript preamble
!  5/20/01  MDG 3.0 f90
! 11/27/01  MDG 3.1 added kind support
! ###################################################################
!
! This file contains a number of Postscript generating subroutines
! which can be used to generate PS-drawings from a Fortran program.
! The routines are based on the Pascal version written by J. Fransaer
! of the Catholic University of Leuven, Belgium.
! Translated into Fortran by M. De Graef on 12/31/92.
!
! The following routines are available (variables are defined in ()s ):
!
!
! ###################################################################
! ###################################################################
!
! A typical program would start as follows
!
!       call PS_openfile(scale)
!       call PS_drawframe(...)
!       ....
!
! and end with 
!
!       call PS_closefile
!
! The code to draw spheres (sp) is taken from the appendix of 
! Earl J. Kirklands book on Advanced Computing in Electron Microscopy,
! and slightly modified to include color.
! ###################################################################

module postscript
!
! All dimensions are in inches 
! pspage      = pagenumber for multi-page files
! psfigwidth  = width of drawing (pagewidth - 2 inches for margins = 6.5)
! psfigheight = height of drawing (pageheight - 2 = 9.0)
! psscale     = scale factor for overall file
! psname      = PostScript file name
! psdash      = used to store a dash pattern
! fonts       = array with standard PostScript fonts (more may be added)
!
! The axis routine has its own common block which has the following entries
! axw         = the width of the square region inside which the entire
!               user space is contained (inches)
! xll,yll     = x and y coordinates (in inches) of the lower left corner of 
!               the square region (measured from lower left corner of page)
!
! The axonometry routine has its own common block which has the following entries

use local

character(55),parameter,private :: PSpreamble(23) = (/ &
        "%!PS-Adobe-3.0                                         ", &
        "%%Creator:                                             ", &
        "%%Title:                                               ", &
        "%%Pages: (atend)                                       ", &
        "%%EndComments                                          ", &
        "/M {moveto} def /N {newpath} def /L {lineto} def       ", &
        "/S {stroke} def /T {translate} def /R {rotate} def     ", &
        "/F {fill} def /Cl {closepath} def                      ", &
        "/circle {N 0 0 1 0 360 arc Cl F} def                   ", &
        "/sp { gsave T scale 1.0 -0.04 0 { pop 3 array astore   ", &
        "{1.02 mul} forall 3 copy setrgbcolor -0.025 0.030 T    ", &
        "circle 0.93 0.93 scale } for                           ", &
        "pop pop pop grestore } def                             ", &
        "/frame {1.0 setgray N left rad add bottom M            ", &
        "right bottom right top rad arcto L right top left top  ", &
        "rad arcto L left top left bottom rad arcto L left      ", &
        "bottom right bottom rad arcto L Cl F 0.0 setgray N     ", &
        "left rad add bottom M right bottom right top rad       ", &
        "arcto L right top left top rad arcto L left top left   ", &
        "bottom rad arcto L left bottom right bottom rad arcto  ", &
        "L Cl S } def                                           ", &
        "%%EndProlog                                            ", &
        "72 dup scale                                           " /)


type postscript_type
 integer(kind=irg)   :: pspage
 real(kind=sgl)      :: psdash(20),psfigwidth,psfigheight,psscale
 character(20)       :: psname
end type

type axonotype
 integer(kind=irg)   :: xi,yi,beta,xmod,ymod,countx,county
 real(kind=sgl)      :: grid,scle,vscle,xstart,ystart
 logical             :: visibility
end type

type axistype
 real(kind=sgl)      :: axw,xll,yll
end type

character(20),parameter :: PSlbl = "Written by MDG, 2001"
character(20),parameter :: PSfonts(5) = (/"Symbol              ", &
                                          "Times-Bold          ", &
                                          "Times-BoldItalic    ", &
                                          "Times-Italic        ", &
                                          "Times-Roman         "/)
type (postscript_type)    :: PS
type (axonotype) :: AXO
type (axistype)  :: AX
integer(kind=irg),allocatable :: imaint(:,:)
integer(kind=irg)             :: imanum

contains 
! ###################################################################
! 
!  subroutine PS_openfile
!
!  Author: Marc De Graef
!  
!  Description: open Postscript file and dump preamble 
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_openfile(dontask)

use local
use io
use files

IMPLICIT NONE

real(kind=sgl)    :: fw, fh
integer(kind=irg) :: i
logical,optional  :: dontask

 imanum = 0
 PS%psfigwidth=6.5; PS%psfigheight=9.0; PS%psscale=1.0
! open file and dump Prolog and Comments sections
 if (present(dontask)) then
! if we get here, it means that we are using a temporary PostScript file
! and we should not go through the regular SafeOpenFile routine.
   open(unit=psunit,file=PS%psname,status='unknown',action='write',form='formatted')
   mess = 'Opening temporary file for PostScript output'; call Message("(A)")
 else
   call SafeOpenFile('ps','formatted',PS%psname)
 end if
 write (psunit,"(A)") PSpreamble(1)
 write (psunit,"(A,' ',A)") trim(PSpreamble(2)), username
 write (psunit,"(A,' ',A)") trim(PSpreamble(3)), progdesc
 do i=4,23
  write (psunit,"(A)") PSpreamble(i)
 end do 
! determine lower left corner and translate to that point
 fw=0.5*(8.50-PS%psscale*PS%psfigwidth)
 fh=0.5*(11.0-PS%psscale*PS%psfigheight)
 write (psunit,"(F12.7,' ',F12.7,' T')") fw,fh
 write (psunit,"(F12.7,' setlinewidth')") 0.01
 write (psunit,"(F12.7,' ',F12.7,' scale')") PS%psscale,PS%psscale
 PS%pspage = 0
end subroutine
! ###################################################################
! 
!  subroutine PS_closefile
!
!  Author: Marc De Graef
!  
!  Description: close Postscript file 
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_closefile

use local
use files

IMPLICIT NONE

 write (psunit,*) 'showpage'
 write (psunit,"(' %%Pages: ',i3)") PS%pspage
 write (psunit,"(' %%EOF')")
 call SafeCloseFile('ps','keep',PS%psname)
end subroutine
! ###################################################################
! 
!  subroutine PS_newpage 
!
!  Author: Marc De Graef
!  
!  Description: start a new page in the PostScript file
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_newpage(frm,btxt)

IMPLICIT NONE

logical        :: frm
character(*)   :: btxt

 if (PS%pspage.ne.0) then
  write (psunit,*) 'showpage saveobj restore'
 end if
 PS%pspage = PS%pspage + 1
 write (psunit,"(' %%Page: ',i3,i3)") PS%pspage-1,PS%pspage
 write (psunit,*) '/saveobj save def'
 call PS_setfont(PSfonts(3),0.18)
 write (psunit,"(1x,F12.7,' ',F12.7,' M (',I8,') show')") 6.75,PS%psfigheight-0.2,PS%pspage
 if (frm.eqv..TRUE.) then
  call PS_drawframe(6.75,PS%psfigheight)
 endif
 call PS_setlinewidth(0.012)
 call PS_textballoon(2.0,9.2,btxt,PSfonts(2),0.25)
 call PS_setfont(PSfonts(5),0.07)
 call PS_text(0.1,-0.1,PSlbl)
end subroutine
! ###################################################################
! 
!  subroutine PS_cellinfo
!
!  Author: Marc De Graef
!  
!  Description: write unit cell information (for drawing programs)
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_cellinfo(xo,yo)

use local
use crystalvars

IMPLICIT NONE

real(kind=sgl)  :: xo,yo

intent(IN)      :: xo,yo

 call PS_setfont(PSfonts(2),0.12/PS%psscale)
 call PS_text(xo,yo,'Filename: '//cell%fname)
 call PS_text(xo,yo-0.16,'a [nm]          ')
 call PS_text(xo,yo-0.30,'b [nm]          ')
 call PS_text(xo,yo-0.44,'c [nm]          ')
 call PS_text(xo,yo-0.58,'alpha [deg]     ')
 call PS_text(xo,yo-0.72,'beta  [deg]     ')
 call PS_text(xo,yo-0.86,'gamma [deg]     ')
 call PS_textvar(xo+0.75,yo-0.16,': ',cell%a)
 call PS_textvar(xo+0.75,yo-0.30,': ',cell%b)
 call PS_textvar(xo+0.75,yo-0.44,': ',cell%c)
 call PS_textvar(xo+0.75,yo-0.58,': ',cell%alpha)
 call PS_textvar(xo+0.75,yo-0.72,': ',cell%beta)
 call PS_textvar(xo+0.75,yo-0.86,': ',cell%gamma)
end subroutine
! ###################################################################
! 
!  subroutine PS_clippath
!
!  Author: Marc De Graef
!  
!  Description: make the last path the clippath
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_clippath

use local

IMPLICIT NONE

 write (psunit,"('Cl clip')")
end subroutine
! ###################################################################
! 
!  subroutine PS_translate
!
!  Author: Marc De Graef
!  
!  Description: redefine the origin of the current coordinate frame
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_translate(x,y)

use local

IMPLICIT NONE

real(kind=sgl)  :: x,y

intent(IN)      :: x,y

 write (psunit,"(F18.7,' ',F18.7,' T')") x,y
end subroutine
! ###################################################################
! 
!  subroutine PS_move    
!
!  Author: Marc De Graef
!  
!  Description: move to a given location
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_move(x,y)

use local

IMPLICIT NONE

real(kind=sgl)  :: x,y

intent(IN)      :: x,y

 write (psunit,"(F18.7,' ',F18.7,' M')") x,y
end subroutine
! ###################################################################
! 
!  subroutine PS_draw 
!
!  Author: Marc De Graef
!  
!  Description: draw a line from the current point to the new point
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_draw(x,y)

use local

IMPLICIT NONE

real(kind=sgl)  :: x,y

intent(IN)      :: x,y

 write (psunit,"(F18.7,' ',F18.7,' L')") x,y
end subroutine
! ###################################################################
! 
!  subroutine PS_setlinewidth
!
!  Author: Marc De Graef
!  
!  Description: set the line width
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_setlinewidth(x)

use local

IMPLICIT NONE

real(kind=sgl)  :: x

intent(IN)      :: x

 write (psunit,"(F12.7,' setlinewidth')") x
end subroutine
! ###################################################################
! 
!  subroutine PS_square
!
!  Author: Marc De Graef
!  
!  Description: draw a square
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_square(x,y,edge)

use local

IMPLICIT NONE

real(kind=sgl)  :: ed,edge,x,y

intent(IN)      :: edge,x,y

 ed=0.5*edge
 write (psunit,"('0.0 setgray')")
 write (psunit,"('newpath')")
 write (psunit,"(2(F12.7,' '),'moveto')") x-ed,y-ed
 write (psunit,"(2(F12.7,' '),'lineto')") x-ed,y+ed
 write (psunit,"(2(F12.7,' '),'lineto')") x+ed,y+ed
 write (psunit,"(2(F12.7,' '),'lineto closepath S')") x+ed,y-ed
end subroutine
! ###################################################################
! 
!  subroutine PS_filledsquare
!
!  Author: Marc De Graef
!  
!  Description: draw a filled square
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_filledsquare(x,y,edge,graylevel)
       
use local

IMPLICIT NONE

real(kind=sgl)  :: ed,edge,x,y,graylevel

intent(IN)      :: edge,x,y,graylevel

 ed=0.5*edge
 write (psunit,"(F12.7,' setgray')") graylevel
 call PS_newpath
 call PS_move(x-ed,y-ed)
 call PS_draw(x-ed,y+ed)
 call PS_draw(x+ed,y+ed)
 write (psunit,"(2(F12.7,' '),'lineto closepath fill S')") x+ed,y-ed
end subroutine
! ###################################################################
! 
!  subroutine PS_cross
!
!  Author: Marc De Graef
!  
!  Description: draw a small cross  
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_cross(x,y,edge,lw)
        
use local

IMPLICIT NONE

real(kind=sgl)  :: ed,edge,x,y,lw

intent(IN)      :: edge,x,y,lw

 ed=0.5*edge
 call PS_setlinewidth(lw)
 call PS_newpath
 call PS_move(x-ed,y-ed)
 call PS_draw(x+ed,y+ed)
 call PS_move(x-ed,y+ed)
 call PS_draw(x+ed,y-ed)
 call PS_stroke
end subroutine
! ###################################################################
! 
!  subroutine PSphere
!
!  Author: Marc De Graef
!  
!  Description: draw a colored sphere
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_sphere(x,y,r,clr)

! method modified from Earl J. Kirkland''s book, page 226, adapted for 
! color PostScript

use local

IMPLICIT NONE

real(kind=sgl)         :: r,x,y
character(3)           :: clr

intent(IN)             :: x,y,r,clr

 if (clr.eq.'red') write (psunit,"(1x,7(f12.5,1x),'sp')") 0.6,0.0,0.0,r,r,x,y
 if (clr.eq.'grn') write (psunit,"(1x,7(f12.5,1x),'sp')") 0.0,0.6,0.0,r,r,x,y
 if (clr.eq.'blu') write (psunit,"(1x,7(f12.5,1x),'sp')") 0.0,0.0,0.6,r,r,x,y
 if (clr.eq.'bro') write (psunit,"(1x,7(f12.5,1x),'sp')") 0.6,0.6,0.4,r,r,x,y
 if (clr.eq.'ylw') write (psunit,"(1x,7(f12.5,1x),'sp')") 0.6,0.6,0.0,r,r,x,y
 if (clr.eq.'pnk') write (psunit,"(1x,7(f12.5,1x),'sp')") 0.6,0.0,0.6,r,r,x,y
 if (clr.eq.'cyn') write (psunit,"(1x,7(f12.5,1x),'sp')") 0.0,0.6,0.6,r,r,x,y
end subroutine
! ###################################################################
! 
!  subroutine PS_arc
!
!  Author: Marc De Graef
!  
!  Description: draw an arc of a circle
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_arc(x0,y0,x,y,radius,ang1,ang2)

use local

IMPLICIT NONE

real(kind=sgl)  :: x0,y0,radius,x,y,ang1,ang2

intent(IN)      :: x0,y0,radius,x,y,ang1,ang2

 write (psunit,"('N ',2(F16.10,' '),' moveto ',5(E16.8,' '),' arc S')") x0,y0,x,y,radius,ang1,ang2
end subroutine
! ###################################################################
! 
!  subroutine PS_circle
!
!  Author: Marc De Graef
!  
!  Description: draw a circle
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_circle(x,y,radius)

use local

IMPLICIT NONE

real(kind=sgl)  :: radius,x,y

intent(IN)      :: radius,x,y

 write (psunit,"('N ',3(F16.10,' '),'0 360 arc Cl S')") x,y,radius
end subroutine
! ###################################################################
! 
!  subroutine PS_filledcircle
!
!  Author: Marc De Graef
!  
!  Description: draw a filled circle
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_filledcircle(x,y,radius,graylevel)
        
use local

IMPLICIT NONE

real(kind=sgl)   :: radius,graylevel,x,y

intent(IN)       :: radius,graylevel,x,y

 write (psunit,"(F12.7,' setgray')") graylevel
 write (psunit,"('N ',3(F12.7,' '),'0 360 arc Cl F')") x,y,radius
end subroutine
! ###################################################################
! 
!  subroutine PS_drawframe
!
!  Author: Marc De Graef
!  
!  Description: draw the main frame
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_drawframe(x,y)
        
use local

IMPLICIT NONE

real(kind=sgl)  :: x,y 

intent(IN)      :: x,y 

 write (psunit,"('N')") 
 call PS_move(0.0,0.0)
 call PS_draw(0.0,y)
 call PS_draw(x,y)
 call PS_draw(x,0.0)
 call PS_draw(0.0,0.0)
 write (psunit,"('Cl S')") 
 write (psunit,"('[0.15 0.03 0.02 0.03] 0 setdash')") 
 write (psunit,"('[] 0 setdash')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_drawrect
!
!  Author: Marc De Graef
!  
!  Description:  draw a rectangle
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_drawrect(x1,y1,x2,y2)
        
use local

IMPLICIT NONE

real(kind=sgl)  :: x1,y1,x2,y2 

intent(IN)      :: x1,y1,x2,y2 

 write (psunit,"('N')") 
 call PS_move(x1,y1)
 call PS_draw(x1,y2)
 call PS_draw(x2,y2)
 call PS_draw(x2,y1)
 call PS_draw(x1,y1)
 write (psunit,"('Cl S')") 
 write (psunit,"('[0.15 0.03 0.02 0.03] 0 setdash')") 
 write (psunit,"('[] 0 setdash')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_line
!
!  Author: Marc De Graef
!  
!  Description: draw a line between two points
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_line(x1,y1,x2,y2)
        
use local

IMPLICIT NONE

real(kind=sgl)  :: x1,y1,x2,y2 

intent(IN)      :: x1,y1,x2,y2 

  call PS_move(x1,y1)
  call PS_draw(x2,y2)
  write (psunit,"('S')")      
end subroutine
! ###################################################################
! 
!  subroutine PS_setdash
!
!  Author: Marc De Graef
!  
!  Description: define a dash pattern
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_setdash(num)

use local

IMPLICIT NONE

integer(kind=irg)  :: i,num

intent(IN)         :: num

 write (psunit,"('[')")
 do i=1,num
  write (psunit,"(F12.7,' ')") PS%psdash(i)
 end do
 write (psunit,"('] ',I4,' setdash')") int(PS%psdash(num+1))
end subroutine
! ###################################################################
! 
!  subroutine PS_closepathS
!
!  Author: Marc De Graef
!  
!  Description: close current path and Stroke
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_closepathS  

use local

IMPLICIT NONE

 write (psunit,"('Cl S')")
end subroutine
! ###################################################################
! 
!  subroutine PS_stroke
!
!  Author: Marc De Graef
!  
!  Description: stroke the current path
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_stroke

use local

IMPLICIT NONE

 write (psunit,"('S ')")
end subroutine
! ###################################################################
! 
!  subroutine PS_gsave
!
!  Author: Marc De Graef
!  
!  Description: save the current graphics settings
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_gsave

use local

IMPLICIT NONE

 write (psunit,"('gsave ')")
end subroutine
! ###################################################################
! 
!  subroutine PS_grestore
!
!  Author: Marc De Graef
!  
!  Description: restore the previous graphics settings
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_grestore

use local

IMPLICIT NONE

 write (psunit,"('grestore ')")
end subroutine
! ###################################################################
! 
!  subroutine PS_closepath
!
!  Author: Marc De Graef
!  
!  Description: close the current path
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_closepath   

use local

IMPLICIT NONE

 write (psunit,"('Cl ')")
end subroutine
! ###################################################################
! 
!  subroutine PS_newpath 
!
!  Author: Marc De Graef
!  
!  Description: start a new path
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_newpath     

use local

IMPLICIT NONE

 write (psunit,"('newpath ')")
end subroutine
! ###################################################################
! 
!  subroutine PS_text
!
!  Author: Marc De Graef
!  
!  Description: draw text at a given location
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_text(x,y,line)

use local

IMPLICIT NONE

real(kind=sgl):: x,y
character(*)  :: line

intent(IN)    :: x,y,line

 write (psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (psunit,"('(',$)") 
 write (psunit,"(A,$)") line
 write (psunit,"(') show')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_textv
!
!  Author: Marc De Graef
!  
!  Description: draw text rotated counterclockwise by 90 degrees
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_textv(x,y,line)

use local

IMPLICIT NONE

real(kind=sgl):: x,y
character(*)  :: line

intent(IN)    :: x,y,line

 write (psunit,"('gsave ')") 
 write (psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (psunit,"('90.0 rotate')") 
 write (psunit,"('( ',A,' ) show')") line
 write (psunit,"('-90.0 rotate grestore')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_texttitle
!
!  Author: Marc De Graef
!  
!  Description: draw the title
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_texttitle(x,y,line,q)

use local

IMPLICIT NONE

real(kind=sgl)     :: x,y,q
character(*)       :: line

intent(IN)         :: x,y,line,q

 write (psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (psunit,"('(',$)") 
 write (psunit,"(A,$)") line
 write (psunit,"('  [x',1PE6.0,'] ) show')") q
end subroutine
! ###################################################################
! 
!  subroutine PS_textvtitle
!
!  Author: Marc De Graef
!  
!  Description: draw a vertical title
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_textvtitle(x,y,line,q)

use local

IMPLICIT NONE

real(kind=sgl)   :: x,y,q
character(*)     :: line

intent(IN)       :: x,y,line,q

 write (psunit,"('gsave ')") 
 write (psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (psunit,"('90.0 rotate')") 
 write (psunit,"('(',$)") 
 write (psunit,"(A,$)") line
 write (psunit,"('  [x',1PE6.0,'] ) show')") q
 write (psunit,"('-90.0 rotate grestore')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_textint 
!
!  Author: Marc De Graef
!  
!  Description: text followed by an integer number
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_textint(x,y,line,val)

use local

IMPLICIT NONE

integer(kind=irg)      :: val
real(kind=sgl)         :: x,y
character(*)           :: line

intent(IN)             :: x,y,line,val

 write (psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (psunit,"('(',$)") 
 write (psunit,"(A,$)") line
 write (psunit,"(I4,$)") val
 write (psunit,"(') show')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_textvar 
!
!  Author: Marc De Graef
!  
!  Description: text followed by a real variable
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
!   5/19/01 MDG 2.0 f90
! ###################################################################
subroutine PS_textvar(x,y,line,val)

use local

IMPLICIT NONE

real(kind=sgl)   :: x,y,val
character(*)     :: line

intent(IN)       :: x,y,line, val

 write (psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (psunit,"('(',$)") 
 write (psunit,"(A,$)") line
 write (psunit,"(F14.4,$)") val
 write (psunit,"(') show')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_textvar8
!
!  Author: Marc De Graef
!  
!  Description: text followed by a double precision variable
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_textvar8(x,y,line,val)

use local

IMPLICIT NONE

real(kind=sgl)   :: x,y
real(kind=dbl)   :: val
character(*)     :: line

intent(IN)       :: x,y,line,val

 write (psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (psunit,"('(',$)") 
 write (psunit,"(A,$)") line
 write (psunit,"(F12.6,$)") val
 write (psunit,"(') show')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_textballoon
!
!  Author: Marc De Graef
!  
!  Description: text inside a rounded balloon
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_textballoon(x,y,line,font,sc)

use local

IMPLICIT NONE

real(kind=sgl):: x,y,sc
character(*)  :: line,font

intent(IN)    :: x,y,line,font,sc

 call PS_setfont(font,sc)
 write (psunit,"('/length (',$)") 
 write (psunit,"(A,$)") line
 write (psunit,"(') stringwidth pop def')") 
 write (psunit,"('/height ',F6.4,' def /border ',F6.4,' def')") 0.11*sc/0.2,0.06*sc/0.2
 write (psunit,"('/bottom ',F12.7,' border sub def')") y
 write (psunit,"('/top ',F12.7,' height add border add def')") y
 write (psunit,"('/left ',F12.7,' border sub def')") x
 write (psunit,"('/right ',F12.7,' length add border add def')") x
 write (psunit,"('/rad 0.04 def frame')") 
 write (psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (psunit,"('(',$)") 
 write (psunit,"(A,$)") line
 write (psunit,"(') show')") 
end subroutine
! ###################################################################
! 
!  subroutine PS_balloon
!
!  Author: Marc De Graef
!  
!  Description: draw a rounded balloon
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_balloon(x,y,le,he,w)

use local

IMPLICIT NONE

real(kind=sgl)  :: x,y,le,he,w

intent(IN)      :: x,y,le,he,w

 write (psunit,"('/he ',F6.4,' def /bo ',F6.4,' def /wi ',F6.4,' def')") he,0.5*w,le
 write (psunit,"('/bottom ',F12.7,' bo add def')") y
 write (psunit,"('/top bottom he add bo sub bo sub def')")  
 write (psunit,"('/left ',F12.7,' bo add def')") x
 write (psunit,"('/right left wi add bo sub def')")  
 write (psunit,"('/rad 0.04 def frame')") 
end subroutine 
! ###################################################################
! 
!  subroutine PS_setfont 
!
!  Author: Marc De Graef
!  
!  Description: select a font and make it active
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_setfont(line,sc)

use local

IMPLICIT NONE

real(kind=sgl):: sc 
character(*)  :: line

intent(IN)    :: line,sc

 write (psunit,"($)") 
 write (psunit,"('/',A,$)") line
 write (psunit,"(' findfont')") 
 write (psunit,"(F6.4,' scalefont ')") sc
 write (psunit,"('setfont')")
end subroutine
! ###################################################################
! 
!  subroutine Printhkl
!
!                                    created: 10/20/98 {9:29:46 AM} 
! 
!  Author: Marc De Graef
!  
!  Description: temporary routine to draw hkl indices in PostScript format
!               used by diffraction program;  to be rewritten to
!               accommodate numbers larger than 10.
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/21/01 MDG 1.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine Printhkl(x,y,h,k,l)

use local

IMPLICIT NONE

character(1),parameter :: numbers(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
integer(kind=irg)      :: h,k,l
real(kind=sgl)         :: x,y,xo,yo,dx,dy,x1,y1
character(1)           :: line

intent(IN)             :: x,y,h,k,l

 call PS_setfont(PSfonts(5),0.065)
 call PS_setlinewidth(0.004)
 xo = 0.050
 yo = 0.050
 dx = 0.050
 dy = 0.065
! THIS ONLY WORKS FOR INDICES -9 <= ... <= 9  !!!
! first index
 x1=x+xo
 y1=y+yo
 line=numbers(abs(h))
 call PS_text(x1,y1,line)
 if (h.lt.0) then
  call PS_line(x1,y1+dy,x1+0.5*dx,y1+dy)
 end if
! second index
 x1=x1+dx
 line=numbers(abs(k))
 call PS_text(x1,y1,line)
 if (k.lt.0) then
  call PS_line(x1,y1+dy,x1+0.5*dx,y1+dy)
 end if
! third index
 x1=x1+dx
 line=numbers(abs(l))
 call PS_text(x1,y1,line)
 if (l.lt.0) then
  call PS_line(x1,y1+dy,x1+0.5*dx,y1+dy)
 end if
end subroutine
! ###################################################################
! 
!  subroutine DumpIndices
!
!                                    created: 10/20/98 {9:29:46 AM} 
! 
!  Author: Marc De Graef
!  
!  Description: draw indices in PostScript format
!               used by stereographic projection programs
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/21/01 MDG 1.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine DumpIndices(S,h,k,l,c,x,y,n)

use local
use crystal

IMPLICIT NONE

character(1),parameter :: numbers(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
logical                :: n
character(1)           :: line,S
integer(kind=irg)      :: c,h,k,l,uvw(3),uvtw(4)
real(kind=sgl)         :: xo,yo,dx,dy,x1,y1,x,y

intent(IN)             :: S,h,k,l,c,x,y,n

 call PS_setfont(PSfonts(5),0.08/PS%psscale)
 xo = 0.050/PS%psscale
 yo = 0.050/PS%psscale
 dx = 0.050/PS%psscale
 dy = 0.075/PS%psscale
 if (n.eqv..FALSE.) then
  xo = -0.30/PS%psscale
  yo = -0.10/PS%psscale
 endif
 if (c.eq.2) then
  xo = -0.30/PS%psscale
  yo = 0.05/PS%psscale
 end if
 uvw =(/ h,k,l /)
 if ((S.eq.'d').AND.(hexset.eqv..TRUE.)) then
  call MilBrav(uvw,uvtw,'34')
  call IndexReduceMB(uvtw)
 end if
! opening bracket
 if (S.eq.'d') then
  call PS_text(x+xo,y+yo,'[')
 else
  call PS_text(x+xo,y+yo,'\(')
        end if
!first index
 x1=x+xo+dx
 y1=y+yo
 if ((S.eq.'d').AND.(hexset.eqv..TRUE.)) then
  line=numbers(abs(uvtw(1)))
 else
  line=numbers(abs(uvw(1)))
 end if
 call PS_text(x1,y1,line)
 if (h.lt.0) then
  call PS_line(x1,y1+dy,x1+0.8*dx,y1+dy)
 end if
! second index
 x1=x1+dx
 if ((S.eq.'d').AND.(hexset.eqv..TRUE.)) then
  line=numbers(abs(uvtw(2)))
 else
  line=numbers(abs(uvw(2)))
 end if
 call PS_text(x1,y1,line)
 if (k.lt.0) then
  call PS_line(x1,y1+dy,x1+0.8*dx,y1+dy)
 end if
! third index (if hexset = .TRUE.) -> put a period
 if (hexset.eqv..TRUE.) then
  x1=x1+dx
  line='.'
  call PS_text(x1,y1,line)
 end if
! last index
 x1=x1+dx
 if ((S.eq.'d').AND.(hexset.eqv..TRUE.)) then
  line=numbers(abs(uvtw(4)))
 else
  line=numbers(abs(uvw(3)))
 end if
 call PS_text(x1,y1,line)
 if (l.lt.0) then
  call PS_line(x1,y1+dy,x1+0.8*dx,y1+dy)
 end if
! closing bracket
 x1=x1+dx
 if (S.eq.'d') then
  call PS_text(x1,y+yo,']')
 else
  call PS_text(x1,y+yo,'\)')
 end if
 dx=dx*0.3
 if (c.eq.2) then
  call PS_line(x+xo,y1-0.02/PS%psscale,x1+dx,y1-0.02/PS%psscale)
 end if
end subroutine
! ###################################################################
! 
!  subroutine PrintIndices
!
!                                    created: 10/20/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: draw indices in PostScript format
!              
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PrintIndices(S,h,k,l,x,y)

use local

IMPLICIT NONE

character(1)     :: S
character(12)    :: line
integer(kind=irg):: h,k,l,hkl(3)
real(kind=sgl)   :: x,y

intent(IN)       :: S,h,k,l,x,y

 hkl = (/ h,k,l /)
 call IndexString(line,hkl,S)
 call PS_text(x,y,line)
end subroutine
! ###################################################################
! 
!  subroutine PS_DumpImage
!
!                                    created: 10/20/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: draw integer image (512x512 maximum size) at given 
!  location with given scale
!              
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_DumpImage(x0,y0,npx,npy,scl)

use local

IMPLICIT NONE

integer(kind=irg)                 :: iq,npx,npy,i,j,ir,iq1,iq2,k
integer(kind=irg),parameter       :: bpp=8
character(2*npx)                  :: bigone
character(3),parameter            :: imnm(20) = (/'i01','i02','i03','i04','i05','i06', &
                                                  'i07','i08','i09','i10','i11','i12','i13','i14', &
                                                  'i15','i16','i17','i18','i19','i20'/)
character(1),parameter            :: hd(0:15) = (/'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'/)
real(kind=sgl)                    :: x0,y0,scl

intent(IN)                        :: x0,y0,npx,npy,scl

 imanum = imanum + 1
! integer image
 call PS_translate(x0,y0)
 write (psunit,"('/picstr ',i3,' string def')") npx
 write (psunit,"(' ')")
 write (psunit,"('/',3A)") (imnm(imanum)(j:j),j=1,3)
 write (psunit,"('{ ',i3,' ',i3,' ',i1,' [',i3,' 0 0 ',i3,' 0 0]')") npx,npy,bpp,npx,npy
 write (psunit,"(' { currentfile picstr readhexstring pop } image } def')")
 write (psunit,"(' gsave ',f7.4,' ',f7.4,' scale ')") scl,scl*float(npy)/float(npx)
 write (psunit,"(3A)") (imnm(imanum)(j:j),j=1,3)
 do j=1,npy
  do i=1,npx
   ir=2*i-1
   iq=imaint(i,j)
   iq1=iq/16
   iq2=mod(iq,16)
   bigone(ir:ir)=hd(iq1)
   ir=ir+1
   bigone(ir:ir)=hd(iq2)
  end do
  k=2*npx
  write (psunit,"(A)") bigone
 end do
 write (psunit,"('grestore')")
 call PS_translate(-x0,-y0)
end subroutine
! ###################################################################
! 
!  subroutine PS_DumpImageDistort
!
!                                    created: 10/20/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: draw integer image (512x512 maximum size) at given 
!  location with given scale which may be different along x and y
!              
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine PS_DumpImageDistort(x0,y0,npx,npy,sclx,scly)

use local

IMPLICIT NONE

integer(kind=irg)                 :: iq,npx,npy,i,j,ir,iq1,iq2,k
integer(kind=irg),parameter       :: bpp=8
character(2*npx)                  :: bigone
character(3),parameter            :: imnm(20) = (/'i01','i02','i03','i04','i05','i06', &
                                                  'i07','i08','i09','i10','i11','i12','i13','i14', &
                                                  'i15','i16','i17','i18','i19','i20'/)
character(1),parameter            :: hd(0:15) = (/'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'/)
real(kind=sgl)                    :: x0,y0,sclx,scly

intent(IN)                        :: x0,y0,npx,npy,sclx,scly

 imanum = imanum + 1
! integer image
 call PS_translate(x0,y0)
 write (psunit,"('/picstr ',i3,' string def')") npx
 write (psunit,"(' ')")
 write (psunit,"('/',3A)") (imnm(imanum)(j:j),j=1,3)
 write (psunit,"('{ ',i3,' ',i3,' ',i1,' [',i3,' 0 0 ',i3,' 0 0]')") npx,npy,bpp,npx,npy
 write (psunit,"(' { currentfile picstr readhexstring pop } image } def')")
 write (psunit,"(' gsave ',f7.4,' ',f7.4,' scale ')") sclx,scly*float(npy)/float(npx)
 write (psunit,"(3A)") (imnm(imanum)(j:j),j=1,3)
 do j=1,npy
  do i=1,npx
   ir=2*i-1
   iq=imaint(i,j)
   iq1=iq/16
   iq2=mod(iq,16)
   bigone(ir:ir)=hd(iq1)
   ir=ir+1
   bigone(ir:ir)=hd(iq2)
  end do
  k=2*npx
  write (psunit,"(A)") bigone
 end do
 write (psunit,"('grestore')")
 call PS_translate(-x0,-y0)
end subroutine
! ###################################################################
! 
!  subroutine IndexReduce
!
!                                    created: 10/20/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: Reduce an index triplet to smallest integers
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine IndexReduce(hkl)

use local

IMPLICIT NONE

integer(kind=irg)  :: hkl(3),mi,i,j
real(kind=sgl)     :: rhkl(3),ir

intent(INOUT)      :: hkl

 mi=100
 do i=1,3
  if ((abs(hkl(i)).lt.mi).and.(hkl(i).ne.0)) mi=abs(hkl(i))
 end do
! then check if this index is a common divider of the others
 j = 0
 do i=1,3
  rhkl(i) = float(hkl(i))/float(mi)
  ir = abs(rhkl(i))-float(int(abs(rhkl(i))))
  if (ir.eq.0.0) j=j+1
 end do
 if (j.eq.3) mi=1
  hkl = int(rhkl*mi)
end subroutine
! ###################################################################
! 
!  subroutine IndexReduceMB
!
!                                    created: 10/20/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: Reduce an index quadruplet to smallest integers
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine IndexReduceMB(hkl)

use local

IMPLICIT NONE

integer(kind=irg)   :: hkl(4),mi,i,j
real(kind=sgl)      :: rhkl(4),ir

intent(INOUT)       :: hkl

 mi=100
 do i=1,4
  if ((abs(hkl(i)).lt.mi).and.(hkl(i).ne.0)) mi=abs(hkl(i))
 end do
! then check if this index is a common divider of the others
 j = 0
 do i=1,4
  rhkl(i) = float(hkl(i))/float(mi)
  ir = abs(rhkl(i))-float(int(abs(rhkl(i))))
  if (ir.eq.0.0) j=j+1
 end do
 if (j.eq.4) mi=1
 hkl = int(rhkl*mi)
end subroutine
! ###################################################################
! 
!  subroutine IndexString
!
!                                    created: 10/21/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: return a string of indices for printing (only deals
!  with indices up to 9)
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/21/98 MDG 1.0 original
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine IndexString(st,hkl,sp)

use local
use crystal 

IMPLICIT NONE

integer(kind=irg)        :: hkl(3),l,hkil(4),i
character(12)            :: st
character(1)             :: sp
character(1),parameter   :: numbers(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)

intent(IN)               :: sp
intent(OUT)              :: st
intent(INOUT)            :: hkl

 do l=1,12
  st(l:l) = ' '
 end do
 l=1
 if (sp.eq.'d') then
  st(l:l)='['
  l=l+1
 else
  st(l:l)='\'
  l=l+1
  st(l:l)='('
  l=l+1
 end if
 if (hexset.eqv..FALSE.) then 
  do i=1,3
   if (hkl(i).lt.0) then
    st(l:l)='-'
   else
    st(l:l)=' '
   end if
   l=l+1
   st(l:l)=numbers(abs(hkl(i)))
   l=l+1
  end do
 else
  if (sp.eq.'d') then 
    call MilBrav(hkl,hkil,'34')
    call IndexReduceMB(hkil)
  else
    hkil(1:2) = hkl(1:2)
    hkil(3) = -(hkl(1)+hkl(2))
    hkil(4) = hkl(3)
  end if
  do i=1,4
   if ((hkl(i).lt.0).and.(i.ne.3)) then
    st(l:l)='-'
   else
    st(l:l)=' '
   end if
   l=l+1
   if (i.eq.3) then 
    st(l:l)='.'
   else
    st(l:l)=numbers(abs(hkil(i)))
   end if
   l=l+1
  end do
 end if
 if (sp.eq.'d') then
  st(l:l)=']'
  l=l+1
  st(l:l)=' '
  l=l+1
  st(l:l)=' '
 else
  st(l:l)='\'
  l=l+1
  st(l:l)=')'
 end if
end subroutine
! ###################################################################
!
!  subroutine DrawSPFrame
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
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine DrawSPFrame(CX,CY,CRad,iview,sp)

use local

IMPLICIT NONE

character(10)     :: instr
character(17)     :: str
real(kind=sgl)    :: CX, CY, CRad 
integer(kind=irg) :: iview(3)
character(1)      :: sp

intent(IN)        :: CX,CY,CRad,sp
intent(INOUT)     :: iview

 call PS_newpage(.FALSE.,'Stereographic Projection')

 call PS_setlinewidth(0.012)
 call PS_circle(CX,CY,CRad)
 call PS_setlinewidth(0.008)
 call PS_line(CX-CRad,CY,CX+CRad,CY)
 call PS_line(CX,CY-CRad,CX,CY+CRad)
 call PS_text(CX-CRad-0.07,CY-0.025,'A')
 call PS_text(CX+CRad+0.03,CY-0.025,'B')
 call PS_text(CX-0.03,CY-CRad-0.08,'M''')
 call PS_text(CX-0.03,CY+CRad+0.05,'M"')
 call PS_cellinfo(0.00,8.30)
 call IndexString(instr,iview,'d')
 call PS_text(CX,8.14,'Viewing Direction '//instr)
 if (sp.eq.'d') then 
  str='direct space'
 else
  str='reciprocal space'
 end if
 call PS_text(CX,8.00,'Projection of '//str)
end subroutine
! ###################################################################
!  subroutine GetIndex
!
!                                    created: 10/13/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: get the u,v,w or h,k,l indices
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/21/01 MDG 2.0 f90
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine GetIndex(ind,sp)

use local
use crystal
use io

IMPLICIT NONE

character(1)      :: sp
integer(kind=irg) :: ind(3),jnd(4),i

intent(IN)        :: sp
intent(OUT)       :: ind

 if (sp.eq.'d') then 
  if (hexset.eqv..FALSE.) then
   mess = 'Enter u,v,w :'; call GetInt(3)
   ind(1:3) = io_int(1:3)
  else
   mess = 'Enter u,v,t,w :'; call GetInt(4)
   jnd(1:4) = io_int(1:4)
   call MilBrav(ind,jnd,'43')
   call IndexReduce(ind)
  end if
 else
  mess = 'Enter h,k,l :'; call GetInt(3)
  ind(1:3) = io_int(1:3)
 endif
end subroutine
end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:  xtalinfo.f90                                                       !
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
!  FILE: "xtalinfo.f90"
!                                    created: 12/11/98 {9:29:46 AM} 
!                                
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: This program generates a multipage PostScript file
!               with a variety of crystallographic output relevant
!               to TEM.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  12/11/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
program xtalinfo

use local
use io
use crystalvars
use crystal
use files
use symmetryvars
use symmetry
use postscript
use diffraction

character(1)     :: sgn,ans
logical          :: first,again,nn,more,topbot
real             :: rr(4),g(3),r(3),cc,M(3,3),negthresh
integer          :: h,k,l,mi,fcnt


 progname = 'xtalinfo.f90'
 progdesc = 'Important crystallographic data for TEM applications'
 call CTEMsoft
 
 inm=2
 SG % SYM_reduce=.TRUE.
 topbot=.FALSE.
! read crystal information
 call CrystalData
 call GetVoltage
 mess = ' Camera Length [mm, R] : '; call GetReal(1)
 camlen = io_real(1)
! generate all atom positions in the fundamental unit cell
 call CalcPositions('v')
! open PostScript file
 call PS_openfile
 pspage = 0
 mess = 'Crystallographic Information'; call Message("(/A/)")
 call InfoPage 
 mess = 'Structure Factors'; call Message("(/A/)")
 call StrucFacPage
 mess = 'Stereographic Projections'; call Message("(/A/)")
 call StereoPage
 mess = 'Diffraction Patterns' ; call Message("(/A/)")
 call DiffPage
!
! close Postscript file
 call PS_closefile
end program
! ###################################################################
!
!  subroutine InfoPage
!
!  Author: Marc De Graef
!
!  Description: creates a page with crystallographic unit cell info
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
subroutine InfoPage

use local
use io
use constants
use crystalvars
use crystal
use files
use symmetryvars
use symmetry
use postscript

logical       :: SaveSpace
real          :: ast,bst,cst,alphast,betast,gammast,coor(5)
integer       :: iap
character(10) :: instr
character(17) :: str
character(1),parameter  :: xsys(7) = (/'c','t','o','h','R','m','a'/)
character(2)  :: brvs
!
! start first page 
 call PS_newpage(.TRUE.,cell % fname)
! write text 
 call PS_text(5.25,0.05,'Distances [nm], angles [degrees]')
 call PS_setlinewidth(0.012)
! direct space information
 call PS_textballoon(0.75,8.45,'Direct Space',PSfonts(2),0.20)
 call PS_setfont(PSfonts(4),0.14)
 call PS_text(1.0,8.00,'a :')
 call PS_text(1.0,7.84,'b :') 
 call PS_text(1.0,7.68,'c :')
 write (psunit,900) 1.25,8.00,cell % a
 write (psunit,900) 1.25,7.84,cell % b
 write (psunit,900) 1.25,7.68,cell % c
 call PS_setfont(PSfonts(1),0.14)
 call PS_text(1.0,7.52,'a :')
 call PS_text(1.0,7.36,'b :')
 call PS_text(1.0,7.20,'g :')
 call PS_setfont(PSfonts(4),0.14)
 write (psunit,900) 1.25,7.52,cell % alpha
 write (psunit,900) 1.25,7.36,cell % beta 
 write (psunit,900) 1.25,7.20,cell % gamma
! space group information and unit cell volume
 call PS_text(3.0,8.00,'Space Group :')
 call PS_setfont(PSfonts(2),0.14)
 write (psunit,905) 4.2,8.00,SYM_SGname(cell % SYM_SGnum),cell % SYM_SGnum
 call PS_setfont(PSfonts(4),0.14)
 call PS_text(3.0,7.84,'Bravais Lattice :')
 call PS_setfont(PSfonts(2),0.14)
 brvs(1:1)=xsys(cell % xtal_system)
 if (cell % xtal_system.ne.5) then 
  brvs(2:2)=SYM_SGname(cell % SYM_SGnum)(2:2)
 else
  brvs(2:2)=' '
 endif
 write (psunit,902) 4.2,7.84,brvs
 call PS_setfont(PSfonts(4),0.14)
 call PS_text(3.0,7.52,'Volume V [nm^3] :')
 write (psunit,900) 4.2,7.52,cell % vol
! reciprocal space information
 call PS_textballoon(0.75,6.75,'Reciprocal Space',PSfonts(2),0.20)
 ast=sqrt(cell % rmt(1,1))
 bst=sqrt(cell % rmt(2,2))
 cst=sqrt(cell % rmt(3,3))
 alphast=acos(cell % rmt(2,3)/bst/cst)*180.0/cPi
 betast =acos(cell % rmt(1,3)/ast/cst)*180.0/cPi
 gammast=acos(cell % rmt(1,2)/ast/bst)*180.0/cPi
 call PS_setfont(PSfonts(4),0.14)
 call PS_text(1.0,6.30,'a* :')
 call PS_text(1.0,6.14,'b* :') 
 call PS_text(1.0,5.98,'c* :')
 write (psunit,900) 1.25,6.30,ast
 write (psunit,900) 1.25,6.14,bst
 write (psunit,900) 1.25,5.98,cst
 call PS_setfont(PSfonts(1),0.14)
 call PS_text(1.0,5.82,'a* :')
 call PS_text(1.0,5.66,'b* :')
 call PS_text(1.0,5.50,'g* :')
 call PS_setfont(PSfonts(4),0.14)
 write (psunit,900) 1.25,5.82,alphast
 write (psunit,900) 1.25,5.66,betast
 write (psunit,900) 1.25,5.50,gammast
! unit cell volume and reciprocal basis vectors
 call PS_text(3.0,6.30,'Volume V* [nm^-3] :')
 write (psunit,900) 4.2,6.30,1.0/cell % vol
! metric tensors and structure matrices
 call PS_textballoon(0.75,5.05,'Important Matrices',PSfonts(2),0.20)
 call DumpMatrix(1.0,4.2,cell % dmt,'g =')
 call DumpMatrix(3.5,4.2,cell % rmt,'g*=')
 call DumpMatrix(1.0,3.2,cell % dsm,'a =')
 call DumpMatrix(3.5,3.2,cell % rsm,'b =')
! asymmetric unit
 call PS_textballoon(0.75,2.80,'Asymmetric Unit',PSfonts(2),0.20)
 dx= 2.75
 dy=-0.16
 topx =-1.75
 topy = 2.3
 x=topx
 y=topy
 db = 0.2
 do i=1,cell % ATOM_ntype
  if (mod(i-1,10).eq.0) then 
   x=x+dx
   call PS_setfont(PSfonts(3),0.14)
   call PS_text(x,y+0.2,'atom ')
   call PS_text(x+0.55,y+0.2,'x ')
   call PS_text(x+0.90,y+0.2,'y ')
   call PS_text(x+1.25,y+0.2,'z ')
   call PS_text(x+1.50,y+0.2,'occ. ')
   call PS_text(x+1.95,y+0.2,'B ')
   call PS_setfont(PSfonts(4),0.12) 
  endif
  y=topy + mod(i-1,10)*dy
  do j=1,3
   coor(j) = cell % ATOM_pos(i,j)
  end do
  coor(4) = cell % ATOM_pos(i,4)
  coor(5) = cell % ATOM_pos(i,5)
  iap = cell % ATOM_type(i)
  call DumpAtom(x,y,ATOM_sym(iap),iap,coor)
 end do
! real number output
 900    format (1x,F12.7,' ',F12.7,' M (',F12.4,') show')
! integer number output
 901    format (1x,F12.7,' ',F12.7,' M (',I8,') show')
! character output
 902    format (1x,F12.7,' ',F12.7,' M (',A,') show')
 905    format (1x,F12.7,' ',F12.7,' M (',A11,'  [#',I3,']) show')
end subroutine
! ###################################################################
!
!  subroutine DumpMatrix
!
!  Author: Marc De Graef
!
!  Description: format a 3x3 matrix in PostScript
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
subroutine DumpMatrix(x,y,g,tt)

use local
use postscript

real(kind=sgl)         :: g(3,3),x,y,dy,dx
character(3) :: tt

! set the matrix label
 call PS_setlinewidth(0.012)
 call PS_setfont(PSfonts(4),0.12) 
 call PS_text(x-0.1,y+0.3,tt)
! draw the left bracket
 dy=0.26
 dx=x+0.1
 call PS_setlinewidth(0.004)
 call PS_line(dx+0.1,y-0.1,dx,y-0.1)
 call PS_line(dx,y-0.1,dx,y+2.5*dy)
 call PS_line(dx,y+2.5*dy,dx+0.1,y+2.5*dy)
 call PS_setlinewidth(0.012)
! write all the entries
 write (psunit,900) dx,y+2.0*dy,g(1,1)
 write (psunit,900) dx,y+1.0*dy,g(2,1)
 write (psunit,900) dx,y+0.0*dy,g(3,1)
 dx=dx+0.5
 write (psunit,900) dx,y+2.0*dy,g(1,2)
 write (psunit,900) dx,y+1.0*dy,g(2,2)
 write (psunit,900) dx,y+0.0*dy,g(3,2)
 dx=dx+0.5
 write (psunit,900) dx,y+2.0*dy,g(1,3)
 write (psunit,900) dx,y+1.0*dy,g(2,3)
 write (psunit,900) dx,y+0.0*dy,g(3,3)
 dx=dx+0.65
! and conclude with the right bracket
 call PS_setlinewidth(0.004)
 call PS_line(dx-0.1,y-0.1,dx,y-0.1)
 call PS_line(dx,y-0.1,dx,y+2.5*dy)
 call PS_line(dx,y+2.5*dy,dx-0.1,y+2.5*dy)
 call PS_setlinewidth(0.012)
!      
 900    format (1x,F12.7,' ',F12.7,' M (',F12.5,') show')
end subroutine
! ###################################################################
!
!  subroutine DumpAtom
!
!  Author: Marc De Graef
!
!  Description: format the atom position and related information
!               for the asymmetric unit.
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
subroutine DumpAtom(x,y,A,AT,coor)

use local
use postscript

real         :: x,y,coor(5)
character(2) :: A
integer      :: AT

 write (psunit,900) x,y,A,AT 
 dx=x+0.40
 write (psunit,910) dx,y,coor(1)
 dx=dx+0.36
 write (psunit,910) dx,y,coor(2)
 dx=dx+0.36
 write (psunit,910) dx,y,coor(3)
 dx=dx+0.36
 write (psunit,910) dx,y,coor(4)
 dx=dx+0.36
 write (psunit,910) dx,y,coor(5)
! formats 
 900    format (1x,F12.7,' ',F12.7,' M (',A2,' [',I2,'] ) show')
 910    format (1x,F12.7,' ',F12.7,' M (',f6.4,') show')
end subroutine
! ###################################################################
!
!  subroutine StrucFacPage
!
!  Author: Marc De Graef
!
!  Description: Structure Factor Information Page
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
subroutine StrucFacPage

use local
use io
use constants
use crystalvars
use crystal
use files
use symmetryvars
use symmetry
use postscript
use diffraction
use dynamical

integer,parameter     :: inm = 4
integer,allocatable   :: family(:,:,:),numfam(:),idx(:)
integer               :: h,k,l,totfam,hh,kk,ll,ind(3),icnt
logical               :: first
logical,allocatable   :: z(:,:,:)
real                  :: g(3),Vmod,Vphase,Vpmod,Vpphase,twopi,vp,vm,vr,vi,xs,qr,qi,bg,rag
real,allocatable      :: Vgg(:,:),ddg(:),gg(:),xi(:),xip(:),twth(:)
character(1)          :: space

! start page 
 call PS_newpage(.TRUE.,'Structure Factor Information')
! write text 
 call PS_text(4.0,0.05,'Distances [nm], angles [mrad], potential [V,rad]')
! label columns
 xs = 0.0
 call PS_setfont(PSfonts(3),0.14)
 call PS_text(xs+1.03,8.15,'h ')
 call PS_text(xs+1.18,8.15,'k ')
 call PS_text(xs+1.33,8.15,'l ')
 call PS_text(xs+1.53,8.15,'p ')
 call PS_text(xs+1.93,8.15,'d ')
 call PS_text(xs+2.38,8.15,'g ')
 call PS_setfont(PSfonts(1),0.14) 
 call PS_text(xs+2.83,8.15,'2q ')
 call PS_setfont(PSfonts(3),0.14) 
 call PS_text(xs+3.33,8.15,'|V| ')
 call PS_text(xs+3.68,8.15,'Vphase ')
 call PS_text(xs+4.18,8.15,'|V''| ')
 call PS_text(xs+4.58,8.15,'V''phase ')
 call PS_setfont(PSfonts(1),0.14) 
 call PS_text(xs+5.23,8.15,'x ')
 call PS_setfont(PSfonts(3),0.14) 
 call PS_text(xs+5.31,8.13,'g ')
 call PS_setfont(PSfonts(1),0.14) 
 call PS_text(xs+5.58,8.15,'x''')
 call PS_setfont(PSfonts(3),0.14) 
 call PS_text(xs+5.66,8.13,'g ')
! initialize parameters
 SG % SYM_reduce=.TRUE.
 thr = 1.E-5
 twopi = 2.0*cPi
 space = 'r'
! allocate all arrays
 allocate(z(-inm:inm,-inm:inm,-inm:inm))
 ii = (2*inm+1)**3
 allocate(family(ii,48,3))
 allocate(numfam(ii))
 allocate(Vgg(ii,6))
 allocate(ddg(ii))
 allocate(gg(ii))
 allocate(xi(ii))
 allocate(xip(ii))
 allocate(twth(ii))
! determine the families of reflections with (hkl)<=(inm,inm,inm)
! first initialize the boolean array z
 z(-inm:inm,-inm:inm,-inm:inm) = .FALSE.
! then loop through all (hkl) values
 first = .TRUE.
 icnt = 1
 totfam=0
 do h=-inm,inm
  ind(1)=-h
  do k=-inm,inm
   ind(2)=-k
   do l=-inm,inm
    ind(3)=-l
! make sure we have not already done this one in another family
    if (.not.z(-h,-k,-l)) then
! if it is a new one, then determine the entire family
     call CalcUcg(ind)
! but ignore the reciprocal lattice point if Vgg is small
     if (rlp%Vmod.ge.thr) then 
! copy family in array and label all its members in z-array
      call CalcFamily(ind,num,space)
      do i=1,num
       do j=1,3
        family(icnt,i,j)=itmp(i,j)
       end do
       z(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
      end do
! store the Fourier coefficient of the lattice potential
      Vgg(icnt,1)=rlp%Vmod
      Vgg(icnt,2)=rlp%Vphase
      Vgg(icnt,3)=rlp%Vpmod
      Vgg(icnt,4)=rlp%Vpphase
      Vgg(icnt,5)=rlp%xgp
      Vgg(icnt,6)=rlp%xg
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
 mess = 'Total number of families        = '; oi_int(1)=icnt; call WriteInt(1,"(I6)")
 mess = 'Total number of family members  = '; oi_int(1)=totfam; call WriteInt(1,"(I6)")
! compute d-spacings, g-spacings, two-theta, and extinction distance
 do k=1,icnt
  do i=1,3
   g(i)=float(family(k,1,i))
  end do
  gg(k)=CalcLength(g,'r')
  ddg(k)=1.0/gg(k)
  twth(k)=2.0*asin(0.5*mLambda*gg(k))
  xip(k) = Vgg(k,5)
  xi(k)  = Vgg(k,6)
 end do
! take care of (0,0,0) reflection
 ind(1)=0
 ind(2)=0
 ind(3)=0
 call CalcUcg(ind)
 Vgg(icnt,1)=rlp%Vmod 
 Vgg(icnt,2)=rlp%Vphase
 Vgg(icnt,3)=rlp%Vpmod 
 Vgg(icnt,4)=rlp%Vpphase
 gg(icnt)=0.0
 twth(icnt)=0.0
 xi(icnt)=0.0
 xip(icnt) = rlp%xgp
! use the spsort.f routine from SLATEC
! to rank by increasing value of gg
 allocate(idx(ii))
 call SPSORT(gg,icnt,idx,1,ier)
! and create output table
 i=1
 j=icnt
 rag = 0.0
 call DumpLine(i,family(j,1,1),family(j,1,2),family(j,1,3),numfam(j),ddg(j),rag,twth(j), &
               Vgg(j,1),Vgg(j,2),Vgg(j,3),Vgg(j,4),xi(j),xip(j))
 do i=2,icnt
  j=idx(i)
  if (ddg(j).ne.0.0) then
   rag = 1.0/ddg(j)
  else
   rag = 0.0
  endif
 call DumpLine(i,family(j,1,1),family(j,1,2),family(j,1,3),numfam(j),ddg(j),rag,twth(j), &
               Vgg(j,1),Vgg(j,2),Vgg(j,3),Vgg(j,4),xi(j),xip(j))
 end do
 deallocate(idx)
 deallocate(z)
 deallocate(family)
 deallocate(numfam)
 deallocate(Vgg)
 deallocate(ddg)
 deallocate(gg)
 deallocate(xi)
 deallocate(xip)
 deallocate(twth)
! integer number output
 901    format (1x,F12.7,' ',F12.7,' M (',I8,') show')
end subroutine
! ###################################################################
!
!  subroutine DumpLine
!
!  Author: Marc De Graef
!
!  Description: formats a line of structure factor information
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
! ###################################################################
subroutine DumpLine(i,h,k,l,p,d,g,th,vm,vp,vpm,vpp,xi,xip)

use local
use postscript

integer         :: h,k,l,p,i
real            :: d,g,th,xi,vm,vp,vpm,vpp,xip

! start a new page (if needed) and label columns
 if ((mod(i-1,50).eq.0).AND.(i.ne.1)) then
  call PS_newpage(.TRUE.,'Structure Factor Information')
! write text 
  call PS_text(4.0,0.05,'Distances [nm], angles [mrad], potential [V,rad]')
! label columns
  xs = 0.0
  call PS_setfont(PSfonts(3),0.14)
  call PS_text(xs+1.03,8.15,'h ')
  call PS_text(xs+1.18,8.15,'k ')
  call PS_text(xs+1.33,8.15,'l ')
  call PS_text(xs+1.53,8.15,'p ')
  call PS_text(xs+1.93,8.15,'d ')
  call PS_text(xs+2.38,8.15,'g ')
  call PS_setfont(PSfonts(1),0.14) 
  call PS_text(xs+2.83,8.15,'2q ')
  call PS_setfont(PSfonts(3),0.14) 
  call PS_text(xs+3.33,8.15,'|V| ')
  call PS_text(xs+3.68,8.15,'Vphase ')
  call PS_text(xs+4.18,8.15,'|V''| ')
  call PS_text(xs+4.58,8.15,'V''phase ')
  call PS_setfont(PSfonts(1),0.14) 
  call PS_text(xs+5.23,8.15,'x ')
  call PS_setfont(PSfonts(3),0.14) 
  call PS_text(xs+5.31,8.13,'g ')
  call PS_setfont(PSfonts(1),0.14) 
  call PS_text(xs+5.58,8.15,'x''')
  call PS_setfont(PSfonts(3),0.14) 
  call PS_text(xs+5.66,8.13,'g ')
 end if
! put all entries in the correct positions
 x = xs+1.0
 y = 8.0 - (mod(i-1,50))*0.15
 call PS_setfont(PSfonts(4),0.12) 
 write (psunit,900) x,y,h       
 x=x+0.15
 write (psunit,900) x,y,k       
 x=x+0.15
 write (psunit,900) x,y,l       
 x=x+0.20
 write (psunit,900) x,y,p       
 x=x+0.25
 write (psunit,901) x,y,d       
 x=x+0.45
 write (psunit,901) x,y,g       
 x=x+0.45
 write (psunit,902) x,y,th*1000.0      
 x=x+0.55
 write (psunit,901) x,y,vm      
 x=x+0.45
 write (psunit,901) x,y,vp      
 x=x+0.45
 write (psunit,901) x,y,vpm     
 x=x+0.45
 write (psunit,901) x,y,vpp     
 x=x+0.45
 if (xi.gt.100000.0) then 
  write (psunit,905) x,y       
 else
  write (psunit,906) x,y,xi       
 endif
 x=x+0.45
 if (xip.gt.100000.0) then 
  write (psunit,905) x,y       
 else
  write (psunit,906) x,y,xip       
 endif
!
 900    format (1x,F12.7,' ',F12.7,' M (',I2,') show')
 901    format (1x,F12.7,' ',F12.7,' M (',F8.5,') show')
 902    format (1x,F12.7,' ',F12.7,' M (',F9.4,') show')
 903    format (1x,F12.7,' ',F12.7,' M (',F12.4,') show')
 904    format (1x,F12.7,' ',F12.7,' M (',I8,') show')
 905    format (1x,F12.7,' ',F12.7,' M ( >100,000 ) show')
 906    format (1x,F12.7,' ',F12.7,' M (',F8.2,') show')
end subroutine
! ###################################################################
!
!  subroutine StereoPage
!
!  Author: Marc De Graef
!
!  Description: creates pages with stereographic projections
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!   4/ 3/01 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
subroutine StereoPage

use local
use symmetryvars
use postscript
use io 
use graphics

character(1)   :: sp
logical        :: topbot
integer        :: hm,km,lm,iview(3)

 SG % SYM_reduce=.TRUE.
 topbot=.TRUE.
 mess = 'Enter the maximum index for h,k and l, or for '; call Message("(A)")
 mess = 'u,v, and w. For a hexagonal system, please use'; call Message("(A)")
 mess = '4-index notation [uv.w] or (hk.l) to determine'; call Message("(A)")
 mess = 'the largest index.'; call Message("(A)")
 mess = 'Enter maximum indices : '; call GetInt(3)
 hm = io_int(1)
 km = io_int(2)
 lm = io_int(3)
! first [001] projection of real space
 sp = 'd'
 iview(1)=0
 iview(2)=0
 iview(3)=1
! call the drawing routine
 call StereoProj(sp,iview,hm,km,lm,topbot)
! then [001] projection of reciprocal space
 sp = 'r'
! call the drawing routine
 call StereoProj(sp,iview,hm,km,lm,topbot)
end subroutine


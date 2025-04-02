! -Fortran-77 Source Code-
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "BWshow.f90"
!                                    created: 4/28/01  {9:29:46 AM} 
!                                last update: 4/28/01 {10:37:32 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: display images and such based on a Bloch wave input file
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  4/28/01  MDG 1.0 original
!  5/27/01  MDG 2.0 f90
! ###################################################################
program BWshow

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript
use error
use dynamical

integer(kind=irg)         :: na,nb,nt,ns,nn
real(kind=sgl)            :: x(2)
character(20)             :: oname
character(2)              :: dtype

 progname = 'BWshow.f90'
 progdesc = 'Display program for Bloch wave data files'
 call CTEMsoft
                                                                               

! get the input filename
 mess = 'Bloch wave input filename : '; call GetStr(oname,20)
! open file, and determine what kind of dataset it is
 open (unit=15,file=oname,form='unformatted',status='old')
 read (15) dtype
 read (15) fname
 select case (dtype)
! Two Beam
  case('TB'); read (15) nn
              read (15) ns
              mess = 'Two Beam data set'; call Message("(A)")
              mess = 'Number of orientations = '; oi_int(1)=ns; call WriteInt(1,"(I3)")
        
! Systematic Row
  case('SR'); read (15) nn
              read (15) ns
              mess = 'Systematic Row data set'; call Message("(A)")
              mess = 'Number of beams        = '; oi_int(1)=nn; call WriteInt(1,"(I3)")
              mess = 'Number of orientations = '; oi_int(1)=ns; call WriteInt(1,"(I3)")
        
! Zone Axis
  case('ZA'); read (15) na,nb
              read (15) ns,nt
              nn = (2*na+1)*(2*nb+1)
              call Message('Zone Axis data set')
  case default; call ErrorMess(8)
 end select
 close (unit=15,status='keep')
! PostScript output file
 call PS_openfile
 PS % pspage=0
!
 if ((dtype.eq.'TB').or.(dtype.eq.'SR')) then
  mess = 'Minimum and maximum foil thickness : '; call GetReal(2)
  mess = 'Number of thicknesses : '; call GetInt(1)
! first, draw grayscale BF and DF images          
  call BWtoPS(nn,ns,io_int(1),io_real(1),io_real(2),oname)
 endif

 call PS_closefile
end program


subroutine BWtoPS(nn,ns,nt,tmin,tmax,oname)

use local 
use io
use constants
use postscript
use graphics

integer(kind=irg)    :: n(1),nn
real(kind=sgl)       :: images(ns,nt,nn),tmin,tmax,y(ns),workx(ns),worky(ns,nn),xmin,xmax,ymin,ymax,imax,imin
integer(kind=irg)    :: ntx,nty,tob
logical              :: npg
character(*)         :: oname
character(12),parameter :: yt(5) = ['gamma^(j)   ', &
                                    'k^(j)_z-k_0 ', &
                                    'alpha^(j)   ', &
                                    'q^(j)       ', &
                                    '            ']
character(40),parameter :: gt(5) = ['Bloch wave eigenvalues                  ', &
                                    'Bloch wave eigenvalues                  ', &
                                    'Bloch wave excitation amplitudes        ', &
                                    'Bloch wave absorption parameters        ', &
                                    'Intensity profiles                      ']
real,parameter          :: xo(4)=[0.5,4.25,0.5,4.25],yo(4)=[5.5,5.5,1.0,1.0]
real,parameter          :: xx(11) = [0.25,2.35,4.45,0.25,2.35,4.45,0.25,2.35,4.45,0.25,2.35], &
                           yy(11) = [6.6,6.6,6.6,4.5,4.5,4.5,2.4,2.4,2.4,0.3,0.3]

 PS % pspage = 0
 call PS_newpage(.FALSE.,'BF/DF Images')
! next compute the beam intensities for a wedge shaped sample
 mess = 'Computing bright field and dark field images'; call Message("(A)")
 call BWtoI(nn,ns,nt,tmin,tmax,oname,images)
! print the images on the top half of the page
 allocate(imaint(ns,nt))
 npx=ns
 npy=nt
 scl = 2.0
 imax=maxval(images)
 imin=minval(images)
 ira = (nn-1)/2
 id  = 6-ira-1
! is it s 2-beam case ?
 if (nn.eq.2) then 
  irmin = 1
  irmax = 2
  id = 0
  lp = 0
! or a multibeam case ?
 else
  irmin = -ira
  irmax = ira
  id = 6-ira-1
  lp = 6
 end if
 do l=irmin,irmax
  j = l+lp
  if ((j.gt.0).and.(j.lt.12)) then
   do i=1,ns
    do ii=1,nt
     imaint(i,ii)=int(255.0*images(i,ii,j-id))
    end do
   end do
   mess = 'dumping image # '; oi_int(1)=j-id; call WriteInt(1,"(I3)")
   call PS_DumpImage(xx(j),yy(j),npx,npy,scl)
  end if
 end do
! next, ask if any other curves are needed
 nf = 0
 isel = 1
 do while (isel.ne.0)
  mess = 'Output of various curves'; call Message("(A)")
  mess = '0. Quit program' ; call Message("(A)")
  mess = '1. Bloch wave eigenvalues (gamma^(j))'; call Message("(A)") 
  mess = '2. Bloch wave eigenvalues (k^(j)_z-k_0)' ; call Message("(A)")
  mess = '3. Bloch wave excitation amplitudes ' ; call Message("(A)")
  mess = '4. Bloch wave absorption parameters ' ; call Message("(A)")
  mess = '5. Store beam intensity vs. thickness' ; call Message("(A)")
  mess = 'Enter your selection :'; call GetInt(1)
  isel = io_int(1)      
  if (isel.ne.0) then
   call ExtractBWdata(workx,worky,kn,isel,ns,nn,oname)
!
   nf = nf+1
   npg = .FALSE.
   if (nf.eq.5) nf=1
   if (nf.eq.1) npg = .TRUE.
   select case (isel) 
    case(1,2,3,4);
     xmin = minval(workx)
     xmax = maxval(workx)
     ymin = minval(worky)
     ymax = maxval(worky)
! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
     if (npg.eqv..TRUE.) then
      call PS_newpage(.FALSE.,'Bloch Wave Results')
     endif
! and call the axis routine to create the plot
! define the drawing location and size
     AX % axw = 4.2
     AX % xll = xo(nf)
     AX % yll = yo(nf)
     y(1:ns) = worky(1:ns,1)
     call axis(ns,workx,y,xmin,xmax,ymin,ymax,.FALSE.,.FALSE.,'lin','lin','CON',1, &
               'BOT','LEF',.FALSE.,.TRUE.,gt(isel),'kt/g',yt(isel))
! superimpose more curves
     do j=2,nn
      y(1:ns) = worky(1:ns,j)
      call axis(ns,workx,y,xmin,xmax,ymin,ymax,.FALSE.,.FALSE.,'lin','lin','CON',1, &
                'BOT','LEF',.TRUE.,.FALSE.,' ',' ',' ')
     end do
    case(5);
        open (dataunit,file='BWshow.out',status='unknown',form='unformatted')
        write (dataunit) images
        close (dataunit,status='keep')
    case default;
   end select
   mess  = 'done'; call Message("(A)")
  end if
 end do
end subroutine
        


subroutine ExtractBWdata(workx,worky,kn,isel,ns,nn,oname)

use local

real(kind=sgl)   :: kttb, kt(3), workx(ns), worky(ns,nn), k(3), kn, kzero
integer(kind=irg):: nn,ns,g(3)
complex(kind=dbl):: W(nn), alph(nn), CG(nn,nn), amp, diag(nn)
character(*)     :: oname
character(2)     :: dtype

 open (unit=15,file=oname,form='unformatted',status = 'old')
 read (15) dtype
 read (15) fname
! Two Beam or Systematic Row
 if ((dtype.eq.'TB').or.(dtype.eq.'SR')) then 
  read (15) i
  read (15) i
  read (15) g
  read (15) k,kzero
  do ik = 1,ns
   read (15) kttb,kn
   read (15) W
   read (15) CG
   read (15) alph
   if (isel.eq.1) then
    workx(ik) = kttb 
    worky(ik,1:nn) = real(W(1:nn))
   end if
   if (isel.eq.2) then
    workx(ik) = kttb 
    worky(ik,1:nn) = -(kn+real(W(1:nn))+kzero)
   end if
   if (isel.eq.3) then
    workx(ik) = kttb 
    worky(ik,1:nn) = abs(alph(1:nn))**2
   end if
   if (isel.eq.4) then
    workx(ik) = kttb 
    worky(ik,1:nn) = imag(W(1:nn))
   end if
  end do
 end if
 close(unit=15,status='keep')
end subroutine


subroutine BWtoI(nn,ns,nt,tmin,tmax,oname,images)

use local
use constants

real(kind=sgl)     :: images(ns,nt,nn),tmin,tmax,kt(3),kttb,kn,kzero,Wr(nn),Wi(nn)
integer(kind=irg)  :: nn,ns,nt,g(3)
complex(kind=dbl)  :: W(nn), alph(nn), CG(nn,nn), amp, diag(nn), q(nn)
character(*)       :: oname
character(15)      :: fname
character(2)       :: dtype

 open (unit=15,file=oname,form='unformatted',status = 'old')
 read (15) dtype
 read (15) fname
! Two Beam
 if ((dtype.eq.'TB').or.(dtype.eq.'SR')) then 
  read (15) i
  read (15) i
  read (15) g
  read (15) kt,kzero
  dz = (tmax-tmin)/float(nt)
  do ik = 1,ns
   read (15) kttb,kn
   read (15) W
   read (15) CG
   read (15) alph
   do i=1,nt
    z = dz*float(i)
    arg = 2.0*sngl(cPi)*z
    Wr = arg*real(W)
    Wi = arg*imag(W)
    q = cmplx(cos(Wr),sin(Wr))
    diag=exp(Wi)*cmplx(cos(Wr),sin(Wr))*alph
    do j=1,nn
     amp = sum(CG(j,1:nn)*diag(1:nn))
     images(ik,i,j) = abs(amp)**2
    end do 
   end do
  end do
 end if
 close (unit=15,status='keep')
end subroutine







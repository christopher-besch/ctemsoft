!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Copyright (C) 2001, 2002    Marc De Graef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    This file is part of CTEMsoft.
!
!    CTEMsoft is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    CTEMsoft is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with CTEMsoft; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "TBBW.f90"
!                                    created: 4/28/01  {9:29:46 AM} 
!                                last update: 4/28/01 {10:37:32 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: Two beam Bloch wave program 
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  4/28/01  MDG 1.0 original
!  5/27/01  MDG 2.0 f90
! ###################################################################
program TBBW

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

integer          :: nt,ns,g(3),n(1),k(3),fn(3)
real             :: ktmax,x(1)
character(20)    :: oname

 progname = 'TBBW.f90'
 progdesc = 'Two-beam bright field-dark field images, using Bloch waves'
 call CTEMsoft
                                                                               

! first get the crystal data and microscope voltage
 SG % SYM_reduce=.TRUE.
 call CrystalData
 call GetVoltage
! generate all atom positions
 call CalcPositions('v')
! get the reciprocal lattice vector
 mess = 'Diffraction vector :'; call Message("(A)")
 call GetIndex(g,'r')
! ask the user for the beam direction 
 mess = 'Beam direction :'; call Message("(A)")
 call GetIndex(k,'d')
 mess = 'Foil normal    :'; call Message("(A)")
 call GetIndex(fn,'d')
! some more parameters
 mess = 'Enter maximum value of k_t in units of g: '; call GetReal(1)
 ktmax = io_real(1)
 mess = 'Number of orientations: '; call GetInt(1)
 ns=io_int(1)
! the file size can be computed as follows:
!  take the number of bytes per variable (4 for real and integer,
!  8 for complex) and count all variables.  Then add 8 bytes for 
!  each write statement (this is an internal Fortran format)
 ifs = (((2+2*2+2)*8+8)+4*8)*ns+((2+15+4+4+3*4+3*4+4)+6*8)
 mess = 'Output file size (bytes) = '; oi_int(1)=ifs; call WriteInt(1,"(I8)")
 mess = 'Output file name : '; call GetStr(oname,20)
! do the computation
 call CalcTBBW(g,float(k),float(fn),ns,ktmax,oname)
end program


subroutine CalcTBBW(g,k,f,ns,ktmax,oname)

use local
use io
use crystal
use crystalvars
use diffraction
use constants
use dynamical

real(kind=sgl)            :: k(3),Vmod,Vphase,Vpmod,Vpphase,pre,upzero,find(3),pr,pi,f(3),kk,l,ktmax,kt(3),kttb,kn,kz,qr,qi,bg
complex(kind=dbl)         :: M(2,2),alph(2),diag(2),amp, CGinv(2,2),Mcp(2,2)
integer(kind=irg)         :: g(3),ind(3),ivec(3),ik,izero, IPIV(2)
character(20)   :: oname

! pre converts from V to U
 pre = 2.0*sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18
! scaling factor for excitation error (2*k_0)
 pre2 = 2.0/sngl(mLambda)
write (*,*) 'pre,pre2 ',pre,pre2
! normal aborption potential Uprime_0
 ind = [0,0,0]
 call CalcUcg(ind)
 Vmod = rlp%Vmod
 Vpmod = rlp%Vpmod
 Vphase = rlp%Vphase
 Vpphase = rlp%Vpphase
 upzero = pre*Vpmod
! tranmitted beam
 izero=1
! determine the dynamical matrix M (all but the diagonal)
! i is the row index
 do i=1,2
  ind = (i-1)*int(g) 
! j is the column index
  do j=1,2
   if (j.ne.i) then
    ivec = ind - (j-1)*int(g)
! use Weickenmeier-Kohl scattering parameters and form factors
    call CalcUcg(ivec)
    M(i,j) = -cPi*cmplx(-aimag(rlp%qg),real(rlp%qg))*cmplx(0.0,1.0/cPi/mLambda)
   end if
  end do
 end do
!
! next we iterate over all incident beam directions, and for
! each direction we complete the M-matrix (diagonal).
!
 dkt = 2.0*ktmax/float(ns-1)
 mess = 'beam tilt step size = '; oi_real(1)=dkt; call WriteReal(1,"(F8.4)")
 find = float(g)
 kk = CalcLength(k,'r')
 gg = CalcLength(find,'r')
 k = k/sngl(mLambda)/kk
 kz = 1.0/mLambda
! open the unformatted output file
 open (unit=15,file=oname,form='unformatted',status ='unknown')
 nn=2
 write (15) 'TB'
 write (15) cell % fname
 write (15) nn
 write (15) ns
 write (15) g
 write (15) k,kz
! loop over the beam directions
 allocate(W(2),CG(2,2))
 do ik = 1,ns
  if (mod(ik,25).eq.0) then
   write (ounit,*) 'completed column ',ik,' of ',ns
  endif
! rescale the wavevector and foil normal
  kt = k + dkt*(float(ik-ns/2)-0.5)*g
  kk = CalcLength(kt,'r')
  kt = kt/sngl(mLambda)/kk
! then complete the diagonal of the M matrix
! i is the row index
  do i=1,nn
   ind = (i-1)*int(g) 
! get the excitation error
   find = float(ind)
   if (i.eq.1) then
    s = 0.0
   else
    s = Calcsg(find,kt,f)
!   write (*,*) s
   endif
! and multiply with 2k_0 and store, along with Uprime_0
   M(i,i) = cmplx(pre2*s,upzero)
  end do
! if (ik.eq.1) then
!  write (*,*) s,kt
!  write (*,*) M
! endif
!
! next, compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
!
! first, make a copy of M, since BWsolve destroys M
  Mcp = M
! then get the eigenvalues and eigenvectors
  call BWsolve(Mcp,W,CG,CGinv,nn,IPIV)
! the alpha coefficients are in the izero column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
  kttb = dkt*(float(ik-ns/2)-0.5)
  kn = -sqrt(kz**2-(kttb*gg)**2)
  W = W/cmplx(2.0*kn,0.0)
  do i=1,nn
   alph(i) = CGinv(i,izero)
  end do
! store eigenvalues in file along with excitation amplitudes, and eigenvector
! matrix (also the wave vector)
  write (15) kttb,kn
  write (15) W
  write (15) CG
  write (15) alph
if (ik.eq.1) then
  write (*,*) kttb,kn
  write (*,*) W
  write (*,*) CG
  write (*,*) alph
endif
 end do
! close the output file
 close(15, status='keep')
 write (ounit,*) 'All data saved in file ',oname
end subroutine

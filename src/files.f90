!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:files.f90                                                            !
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
!  FILE: "files.f90"
!                                    created: 1/5/99 {11:26:49 AM} 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: everything that has to do with file-based input-output
! 
!  History 
!
! 
!  modified by  rev reason
!  -------- --- --- -----------
!    1/5/99 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
!

module files

contains
! ###################################################################
! 
!  subroutine ResetData
!
!                                    created: 10/13/98 {9:29:46 AM} 
!
!  Author: Marc De Graef
!  
!  Description: reset all unitcell and symmetry variables to zero
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine ResetData

use local
use crystalvars
use symmetryvars

IMPLICIT NONE

! initialize cell 
 cell%a = 0.0_sgl
 cell%b = 0.0_sgl
 cell%c = 0.0_sgl
 cell%alpha = 0.0_sgl
 cell%beta  = 0.0_sgl
 cell%gamma = 0.0_sgl
 cell%vol   = 0.0_sgl
 cell%dmt = 0.0_dbl
 cell%rmt = 0.0_dbl
 cell%dsm = 0.0_dbl
 cell%rsm = 0.0_dbl
 cell%krdel = 0.0_dbl
 cell%ATOM_type = 0
 cell%ATOM_ntype = 0
 cell%SYM_SGnum = 0
 cell%xtal_system = 0
 cell%SYM_SGset = 0
 cell%ATOM_pos = 0.0_dbl
 cell%fname = '               '
 cell%SYM_reduce = .FALSE.
 cell%SYM_second = .FALSE.
 cell%SYM_trigonal = .FALSE.

! initialize all symmetry variables
 SG%SYM_GENnum = 0
 SG%SYM_MATnum = 0
 SG%SYM_NUMpt  = 0
 SG%SYM_reduce = .FALSE.
 SG%SYM_trigonal = .FALSE.
 SG%SYM_second = .FALSE.
 SG%SYM_centrosym = .FALSE. 
 SG%SYM_c = 0.0_dbl
 SG%SYM_data = 0.0_dbl
 SG%SYM_direc = 0.0_dbl
 SG%SYM_recip = 0.0_dbl
 SG%SYM_name = '           '
end subroutine

! ###################################################################
! 
!  subroutine DumpXtalInfo
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: Write a brief summary of the crystal structure on the screen
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine DumpXtalInfo    

use local
use io
use crystalvars
use symmetryvars

IMPLICIT NONE

 mess = 'Crystal Structure Information'; call Message("('-->',A,'<--')")
 mess = 'a [nm]             : '; oi_real(1) = cell%a; call WriteReal(1,"(2x,F8.4)")
 mess = 'b [nm]             : '; oi_real(1) = cell%b; call WriteReal(1,"(2x,F8.4)")
 mess = 'c [nm]             : '; oi_real(1) = cell%c; call WriteReal(1,"(2x,F8.4)")
 mess = 'alpha [deg]        : '; oi_real(1) = cell%alpha; call WriteReal(1,"(2x,F8.4)")
 mess = 'beta  [deg]        : '; oi_real(1) = cell%beta ; call WriteReal(1,"(2x,F8.4)")
 mess = 'gamma [deg]        : '; oi_real(1) = cell%gamma; call WriteReal(1,"(2x,F8.4)")
 mess = 'Space group        : #'; oi_int(1) = cell%SYM_SGnum; call WriteInt(1,"(2x,I4)")
 mess = 'Space group symbol : ';  call WriteStr(SYM_SGname(cell%SYM_SGnum))
 mess = 'Generator String   : ';  call WriteStr(SYM_GL(cell%SYM_SGnum))
 if ((cell%SYM_SGset.eq.2).AND.(cell%xtal_system.ne.5)) then 
  mess = 'Using second origin setting'; call Message("(A)")
 endif
 if ((cell%SYM_SGset.eq.2).AND.(cell%xtal_system.eq.5)) then 
  mess = 'Using rhombohedral parameters'; call Message("(A)")
 endif
end subroutine
!
! ###################################################################
! 
!  subroutine CrystalData
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 10/20/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: load or generate crystal data
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine CrystalData(fname)

use local
use io
use crystalvars
use crystal
use symmetry

IMPLICIT NONE

character(15),OPTIONAL,INTENT(IN)  :: fname

 if (PRESENT(fname)) then 
  strucdef = .TRUE.
  cell%fname = fname
 else
  strucdef = .FALSE.
  mess = 'Load file (0) or new data (1) ? '; call GetInt(1)
 end if
 if (.not.strucdef) then
  if (io_int(1).ne.0) then
   cell%SYM_SGset=0
   call GetLatParm
   call CalcMatrices
   call GetSpaceGroup
   call GenerateSymmetry(.TRUE.)
   call GetAsymPos
   call SaveData
   call ResetData
  end if
 end if
 call ReadDataFile
end subroutine 
!
! ###################################################################
! 
! function SaveData
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 12/14/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: save crystal structure data to file
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine SaveData

use local
use io
use crystalvars
use crystal
 
IMPLICIT NONE

 call SafeOpenFile('xt','unformatted',cell%fname)
 open (dataunit,file=cell%fname,status='unknown',form='unformatted')
! save lattice parameters, crystal system, space group number and contents
! of the asymmetric unit.
 write (dataunit) cell%xtal_system, cell%a,cell%b,cell%c,cell%alpha,cell%beta,cell%gamma
 write (dataunit) cell%ATOM_type, cell%ATOM_pos, cell%ATOM_ntype, cell%SYM_SGnum, cell%SYM_SGset
 call SafeCloseFile('xt','keep',cell%fname)
 mess = 'Crystal data saved in file '//cell%fname; call Message("(A)")
 strucdef = .FALSE.
end subroutine
!
! ###################################################################
! 
! subroutine ReadDataFile
!
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 12/14/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: load crystal data in memory
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine ReadDataFile(fr)

use local
use io
use crystalvars
use crystal
use symmetry
use symmetryvars

IMPLICIT NONE

integer(kind=irg) :: i,ipg,isave
logical           :: fread = .TRUE.
logical,optional  :: fr

 if (present(fr)) fread=.FALSE.
 call SafeOpenFile('xt','unformatted',cell%fname,fread)
 read (dataunit) cell%xtal_system, cell%a,cell%b,cell%c,cell%alpha,cell%beta,cell%gamma
 read (dataunit) cell%ATOM_type, cell%ATOM_pos, cell%ATOM_ntype, cell%SYM_SGnum, cell%SYM_SGset
 call SafeCloseFile('xt','keep',cell%fname)
 strucdef = .TRUE.
 hexset = .FALSE.
 if (cell%xtal_system.eq.4) hexset = .TRUE.
 if ((cell%xtal_system.eq.5).AND.(cell%SYM_SGset.ne.2)) hexset = .TRUE.
 call DumpXtalInfo
 call CalcMatrices
! First generate the point symmetry matrices, then the actual space group.
! Get the symmorphic space group corresponding to the point group
! of the actual space group
 ipg=0
 do i=1,32
  if (SGPG(i).le.cell%SYM_SGnum) ipg=i
 end do
! if the actual group is also the symmorphic group, then both 
! steps can be done simultaneously, otherwise two calls to 
! GenerateSymmetry are needed.
 if (SGPG(ipg).eq.cell%SYM_SGnum) then
  call GenerateSymmetry(.TRUE.)
 else
  isave = cell%SYM_SGnum
  cell%SYM_SGnum = SGPG(ipg)
  call GenerateSymmetry(.TRUE.)
  cell%SYM_SGnum = isave
  call GenerateSymmetry(.FALSE.)
 end if
end subroutine
! ###################################################################
! 
! function SafeOpenFile
! 
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 12/14/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: open a data file, based on example 10-1 in 
!  Stephen Chapman''s Fortran 90/95 book.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine SafeOpenFile(ftyp,frm,fname,fread)

use local
use io


IMPLICIT NONE

character(1)        :: ans
character(2)        :: ftyp
character(20)       :: fname
character(*)        :: frm
logical,OPTIONAL,INTENT(IN) :: fread
logical             :: lexist,lopen 
integer             :: iunit

intent(IN)          :: ftyp,frm
intent(OUT)         :: fname

 select case (ftyp)
  case('xt'); iunit = dataunit
  case('d1'); iunit = dataunit
  case('d2'); iunit = dataunit2
  case('d3'); iunit = dataunit3
  case('ps'); iunit = psunit
 end select
 lopen = .FALSE.
 openfile: do while(.not.lopen)
! get filename
  if (PRESENT(fread)) then
   if (fread.eqv..TRUE.) then 
    mess = ' Enter input file name : '; call GetStr(fname,20)
   end if
  else
   mess = ' Enter output file name : '; call GetStr(fname,20)
  end if
! does file already exist ?
  inquire(file=fname,exist=lexist)
  exists: if (.not.lexist) then
! ok, file does not already exist, so open it
    open(unit=iunit,file=fname,status='new',action='write',form=frm)
    lopen = .TRUE.
  else
! file exists.  Should we replace it or read it ?
   if (PRESENT(fread)) then
     open(unit=iunit,file=fname,status='old',action='read',form = frm)
     lopen = .TRUE.
   else
    mess = ' Output file exists.  Overwrite it ? (y/n) '; call GetStr(ans,1)
    replace: if ((ans == 'Y').or.(ans == 'y')) then
! open file
     open(unit=iunit,file=fname,status='replace',action='write',form = frm)
     lopen = .TRUE.
    end if replace
   end if
  end if exists
 end do openfile
 if (lopen) then
  if (PRESENT(fread)) then
   mess = ' File '//fname//' open for input '; call Message("(A)")
  else
   mess = ' File '//fname//' open for output '; call Message("(A)")
  end if
 end if
end subroutine
! ###################################################################
! 
! function SafeCloseFile
! 
!                                    created: 10/20/98 {9:29:46 AM} 
!                                last update: 12/14/98 {9:29:46 AM} 
!  Author: Marc De Graef
!  
!  Description: close data file
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/20/98 MDG 1.0 original
!   5/19/01 MDG 2.0 f90 version
!  11/27/01 MDG 2.1 added kind support
! ###################################################################
subroutine SafeCloseFile(ftyp,stt,fname,quiet)

use local
use io
use error

IMPLICIT NONE

integer(kind=irg)   :: iunit,ier
character(2)        :: ftyp
character(*)        :: stt
character(20)       :: fname
logical,OPTIONAL,INTENT(IN) :: quiet

intent(IN)          :: ftyp,stt,fname


 select case (ftyp)
  case('xt'); iunit = dataunit
  case('d1'); iunit = dataunit
  case('d2'); iunit = dataunit2
  case('d3'); iunit = dataunit3
  case('ps'); iunit = psunit
 end select
 close(unit=iunit,status=stt,iostat=ier)
 if (ier.ne.0) call FatalError('SafeCloseFile: data file NOT saved',' ')
 if (.not.PRESENT(quiet)) then
  mess = ' Data has been stored in file '//fname; call Message("(A)")
 end if
end subroutine


end module

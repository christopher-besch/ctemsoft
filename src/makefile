# ###################################################################
#  makefile     
#
#                                    created:  5/19/01 
#  Author: Marc De Graef
#  
#  Description:  General makefile for unix version of CTEMsoft package.
#
#                This limited version contains the main CTEM library
#                and the program files for Chapters 1 through 4.
# 
#                Typing    make all   should do the trick.  
#                Use       make clean  to clean up all files
#
#  History
# 
#  modified by  rev reason
#  -------- --- --- -----------
#  11/05/03 MDG 1.0 original
# ###################################################################
#
##########################################################################
# USER CONFIGURABLE PARAMETERS # 
################################
# path to the rm command
RM = /bin/rm -f
##########################################################################
# NO USER CONFIGURABLE PARAMETERS BELOW #
#########################################
LIBTEM = libtem.a
CHAP1= chap1
CHAP2= chap2
CHAP3= chap3
CHAP4= chap4
MAKE = make -f m.routines
SHELL = /bin/sh
EXE = ../exe
#
# create the libtem.a library
#
F90FILES = local.f90 io.f90 error.f90 constants.f90 math.f90 crystalvars.f90 symmetryvars.f90 crystal.f90 \
         symmetry.f90 files.f90 others.f90 postscript.f90 diffraction.f90 tiff.f90 ctf.f90 \
         lorentz.f90 graphics.f90 timing.f90 rksuite_90.f90
F77FILES = eispack.f spsort.f

O90FILES = $(F90FILES:.f90=.o)
O77FILES = $(F77FILES:.f=.o)

$(LIBTEM): 
	for i in $(O77FILES); do \
	  ($(MAKE) $$i ); \
	done
	for i in $(O90FILES); do \
	  ($(MAKE) $$i); \
	done
	ar -vr $(LIBTEM) $(O90FILES) $(O77FILES)

#
# general copyright statement program
#
GPL: 
	$(MAKE) GPL 

# here are all the program names, ranked by chapter
F90_C1_FILES = star orbit listSG latgeom family stereo drawcell orient    
# F90_C2_FILES = vr vrfft
F90_C2_FILES = vr
F90_C3_FILES = lens zap holz    
F90_C4_FILES = xtalinfo
F90_C6_FILES = BWshow qg TBBFDF TBBW TBSM

# Chapter 1
$(CHAP1): 
	for i in $(F90_C1_FILES); do \
	  ($(MAKE) $$i); \
	done

# Chapter 2
$(CHAP2): 
	for i in $(F90_C2_FILES); do \
	  ($(MAKE) $$i); \
	done

# Chapter 3
$(CHAP3): 
	for i in $(F90_C3_FILES); do \
	  ($(MAKE) $$i); \
	done

# Chapter 4
$(CHAP4): 
	for i in $(F90_C4_FILES); do \
	  ($(MAKE) $$i); \
	done

# Chapter 6
$(CHAP6):
	for i in $(F90_C6_FILES); do \
          ($(MAKE) $$i); \
        done

all: $(LIBTEM) GPL $(CHAP1) $(CHAP2) $(CHAP3) $(CHAP4) $(CHAP6)

clean	: 
	@echo "removing all executables"
	$(RM) $(EXE)/*
	@echo "removing all mod files"
	$(RM) *.mod
	@echo "removing all object files"
	$(RM) *.o
	@echo "removing all libraries"
	$(RM) *.a



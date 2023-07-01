HOST = gcc
#HOST = gcc-openmp

PROJECT1 = int2-laprouts3d
PROJECT2 = int2-test-lfmm3d-mps

# FC - fortran compiler
# FFLAGS - fortran compiler flags

ifeq ($(HOST),gcc)
    FC=gfortran -w
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy
endif

# Test objects
#
COM = ../../src/Common
LAP = ../../src/Laplace

.PHONY: all clean

default: all laprouts3d.c mexfile


OBJECTS = lfmm3d_mps0.o \
    pts_tree3d0.o \
    $(COM)/hkrand.o \
    $(COM)/dlaran.o \
    $(COM)/prini.o \
    $(COM)/rotgen.o \
    $(COM)/legeexps.o \
    $(COM)/rotviarecur.o \
    $(COM)/yrecursion.o \
    $(COM)/rotproj.o \
    $(COM)/dfft.o \
    $(COM)/fmmcommon.o \
    $(COM)/tree_routs3d.o \
    $(COM)/cumsum.o \
    $(COM)/pts_tree3d.o \
    $(LAP)/l3dterms.o \
    $(LAP)/l3dtrans.o \
    $(LAP)/laprouts3d.o \
    $(LAP)/lapkernels.o \
    $(LAP)/lpwrouts.o \
    $(LAP)/lndiv.o \
    $(LAP)/lwtsexp_sep1.o \
    $(LAP)/lfmm3d.o \

# mex Laplace routine for debugging
MWRAP = ~/mwrap/mwrap
FDIR=$$(dirname `gfortran --print-file-name libgfortran.dylib`)
MFLAGS += -L${FDIR}
MEX = $(shell ls -d /Applications/MATLAB_R20**.app)/bin/mex
laprouts3d.c: laprouts3d.mw
	$(MWRAP) -c99complex -mex laprouts3d -mb -list laprouts3d.mw
	$(MWRAP) -c99complex -mex laprouts3d -c laprouts3d.c laprouts3d.mw
mexfile: laprouts3d.c lfmm3d_mps0.o pts_tree3d0.o $(LAP)/lpwrouts.o $(LAP)/lndiv.o $(LAP)/lwtsexp_sep1.o $(LAP)/lapkernels.o $(LAP)/laprouts3d.o $(LAP)/l3dterms.o $(LAP)/l3dtrans.o $(COM)/fmmcommon.o $(COM)/prini.o $(COM)/yrecursion.o $(COM)/rotproj.o $(COM)/rotviarecur.o $(COM)/dfft.o $(COM)/hkrand.o $(COM)/dlaran.o $(COM)/tree_routs3d.o $(COM)/cumsum.o $(COM)/pts_tree3d.o $(COM)/rotgen.o
	$(MEX) laprouts3d.c lfmm3d_mps0.o pts_tree3d0.o $(LAP)/lpwrouts.o $(LAP)/lndiv.o $(LAP)/lwtsexp_sep1.o $(LAP)/lapkernels.o $(LAP)/laprouts3d.o $(LAP)/l3dterms.o $(LAP)/l3dtrans.o $(COM)/fmmcommon.o $(COM)/prini.o $(COM)/yrecursion.o $(COM)/rotproj.o $(COM)/rotviarecur.o $(COM)/dfft.o $(COM)/hkrand.o $(COM)/dlaran.o $(COM)/tree_routs3d.o $(COM)/cumsum.o $(COM)/pts_tree3d.o $(COM)/rotgen.o -largeArrayDims $(MFLAGS) -lgfortran -lm -lstdc++
#mexfile: laprouts3d.c lfmm3d_mps0.o pts_tree3d0.o $(LAP)/lfmm3d.o $(LAP)/lpwrouts.o $(LAP)/lndiv.o $(LAP)/lwtsexp_sep1.o $(LAP)/lapkernels.o $(LAP)/laprouts3d.o $(LAP)/l3dterms.o $(LAP)/l3dtrans.o $(COM)/fmmcommon.o $(COM)/prini.o $(COM)/yrecursion.o $(COM)/rotproj.o $(COM)/rotviarecur.o $(COM)/dfft.o $(COM)/hkrand.o $(COM)/dlaran.o $(COM)/tree_routs3d.o $(COM)/cumsum.o $(COM)/pts_tree3d.o $(COM)/rotgen.o
#	$(MEX) laprouts3d.c lfmm3d_mps0.o pts_tree3d0.o $(LAP)/lfmm3d.o $(LAP)/lpwrouts.o $(LAP)/lndiv.o $(LAP)/lwtsexp_sep1.o $(LAP)/lapkernels.o $(LAP)/laprouts3d.o $(LAP)/l3dterms.o $(LAP)/l3dtrans.o $(COM)/fmmcommon.o $(COM)/prini.o $(COM)/yrecursion.o $(COM)/rotproj.o $(COM)/rotviarecur.o $(COM)/dfft.o $(COM)/hkrand.o $(COM)/dlaran.o $(COM)/tree_routs3d.o $(COM)/cumsum.o $(COM)/pts_tree3d.o $(COM)/rotgen.o -largeArrayDims $(MFLAGS) -lgfortran -lm -lgomp -lstdc++

all: $(OBJECTS) 
	$(FC) $(FFLAGS) test_laprouts3d.f -o $(PROJECT1) $(OBJECTS) 
	$(FC) $(FFLAGS) test_lfmm3d_mps.f -o $(PROJECT2) $(OBJECTS) 
	./$(PROJECT1)
	./$(PROJECT2)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

clean: 
	rm -f $(OBJECTS) $(PROJECT1) $(PROJECT2) fort.13

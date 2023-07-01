#
#  Generate Matlab/Octave MEX files
#
#  Date: November 9, 2017
#

###TARGET=matlab-maci64
TARGET=matlab-maci64-openmp

##############################################################################
#
#  Matlab | MacOS 64-bit | gfortran compiler
#
##############################################################################


ifeq ($(TARGET), matlab-maci64)

FFLAGS=-g -O2 -fPIC -fallow-argument-mismatch -std=legacy

# Use the following for Matlab + gfortran
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
MEX=   $(shell ls -d /Applications/MATLAB_R202*.app)/bin/mex $(OOFLAG) CC=gcc LD=gcc ARCH=maci64
# Use the following for 64-bit MEX
###MEX=   mex -largeArrayDims $(OOFLAG) ARCH=maci64
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap/mwrap
MWFLAGS=-c99complex -mex        
MEX_LIBS=-lgfortran -lm
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), matlab-maci64-openmp)

FFLAGS=-g -O2 -fPIC --openmp -fallow-argument-mismatch -std=legacy

# Use the following for Matlab + gfortran
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
MEX=   $(shell ls -d /Applications/MATLAB_R202*.app)/bin/mex $(OOFLAG) CC=gcc LD=gcc ARCH=maci64
# Use the following for 64-bit MEX
###MEX=   mex -largeArrayDims $(OOFLAG) ARCH=maci64
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap/mwrap
MWFLAGS=-c99complex -mex        
MEX_LIBS=-lgfortran -lgomp -lm
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif

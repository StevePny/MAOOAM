
ifndef COMPILER
  COMPILER=gfortran
endif
ifdef INTMODE
  # Define integrator module
  INTEGRATOR = $(INTMODE)_integrator_decoupled
  TL_AD_INTEGRATOR = $(INTMODE)_tl_ad_integrator_decoupled
else
  # Define integrator module
  INTEGRATOR = rk2_integrator_decoupled
  TL_AD_INTEGRATOR = rk2_tl_ad_integrator_decoupled
endif
ifdef RES
  ST=_$(RES)
  $(warning Code will $(RES) the inner products in arrays! It may induce a huge memory usage!)
endif
# Define relevant compiler options for debug case and normal case for each compiler
# Normal Case:
IFORT_COMPILER_FLAGS = -assume byterecl -O2
IFORT_LD_FLAGS = -mkl
GCC_COMPILER_FLAGS = -O2 -Wall
GCC_LD_FLAGS = -llapack -lblas
# Debug Case:
IFORT_COMPILER_FLAGS_DBG = -O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl
IFORT_LD_FLAGS_DBG = -llapack -lblas
GCC_COMPILER_FLAGS_DBG = -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow 
GCC_LD_FLAGS_DBG = -llapack -lblas


SHELL = /bin/sh
DEBUG = false
#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif


ifeq ($(COMPILER),gfortran)
  #Derive gfortran version from gcc version and check if >=4.6.0.
  GVERSION := $(shell gcc -dumpversion)
  GE40600 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40600)
  
  ifeq ($(COMPILER),gfortran)
    ifeq ($(GE40600), 0)
#     $(error GCC/gfortran version $(GVERSION) detected. This program uses features of gfortran >= 4.6.0. Please update your toolchain)
    endif
  else
    $(warning This code was optimized for gfortran >= 4.6.0 and has not been tested with other compilers)
  endif
  COMPILER_FLAGS = $(GCC_COMPILER_FLAGS)
  LD_FLAGS = $(GCC_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(GCC_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(GCC_LD_FLAGS_DBG)
  IFPORT = ifport.mod
else ifeq ($(COMPILER),ifort)
  $(warning This code was optimized for ifort >= 14.0.2 and has not been tested with other compilers)
  COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
  LD_FLAGS = $(IFORT_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(IFORT_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(IFORT_LD_FLAGS_DBG)
else
  $(warning This code was optimized for gfortran >= 4.6.0 and ifort >= 14.0.2 and has not been tested with other compilers)
endif

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

MODULES = util.mod params.mod inprod_analytic$(ST).mod aotensor_def.mod ic_def.mod $(INTEGRATOR).mod stat.mod tensor.mod tl_ad_tensor.mod $(TL_AD_INTEGRATOR).mod $(IFPORT) maooam_atmos_wrapper.mod 
MODULE_OBJECTS = $(MODULES:.mod=.o)

all: maooam_atmos_wrapper.o archive

debug: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
debug: LD_FLAGS = $(LD_FLAGS_DBG)
debug: all

%.mod: %.f90
	$(COMPILER) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%.o: %.f90 $(MODULES)
	$(COMPILER) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%: %.o $(MODULE_OBJECTS)
	$(COMPILER) $(COMPILER_FLAGS) $^ -o $@ $(LD_FLAGS) 

inprod_analytic$(ST).mod: params.mod util.mod
aotensor_def.mod: tensor.mod inprod_analytic$(ST).mod
$(INTEGRATOR).mod: aotensor_def.mod
ic_def.mod: util.mod aotensor_def.mod
stat.mod: params.mod
tl_ad_tensor.mod: aotensor_def.mod
$(TL_AD_INTEGRATOR).mod: tl_ad_tensor.mod
util.mod: $(IFPORT)
#STEVE: adding to support compilation of wrapper module
maooam_wrapper.mod: $(MODULES)

#STEVE:
archive:
	ar cr libmaooam_atmos.a *.o
	rm *.o

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test

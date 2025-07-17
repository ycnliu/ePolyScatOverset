#
#
#!/bin/bash
ifndef pe
pe := $(CURDIR)
endif
ifndef pt
pt := $(CURDIR)/tests
export pt
endif
ifndef od
od := $(pe)/tests/outdir.$(MACH)
export od
endif
ifndef po
po := $(pe)/tests/outdir.$(MACH)
export po
endif
# Compiler selection (edit as needed)
FORTMPI  := mpif90
CXX      := mpicxx
NVCC     := nvcc

FFLAGS   := -g -O2 -mtune=native -fbounds-check -fallow-argument-mismatch -fcheck=all
CXXFLAGS := -g -O2 -std=c++17 -march=native
NVCCFLAGS := -O2 -arch=sm_60

# Linker flags (edit for your environment)
LFLAGS := -lcudart -lstdc++ 

IGNORE := $(shell mkdir $(od) >& /dev/null )
#
ifndef COMPILER
machcomp := $(MACH)
else
machcomp := $(MACH)$(COMPILER)
endif
bindir := $(pe)/bin/$(machcomp)
objdir := $(pe)/obj/$(machcomp)
binrel := bin/$(machcomp)
objrel := obj/$(machcomp)
IGNORE := $(shell mkdir $(bindir) >& /dev/null )
IGNORE := $(shell mkdir $(objdir) >& /dev/null )

srcdir := $(pe)/src
GITSHA1 := $(shell git log -1 --format=%H)

export pe bindir
GITDATE := $(shell git log -1 --format=%cD)
export machcomp
export bindir
export objdir
export srcdir
export GITSHA1
export GITDATE

TESTS := $(basename $(notdir $(wildcard tests/test[0-9][0-9].inp)))
TESTSCKOHN := $(basename $(notdir $(wildcard tests/test[5-9][0-9].inp)))
TESTSX := $(basename $(notdir $(wildcard tests/testx*.inp)))

#
#
include include/$(machcomp).mk
#
# use default target for running tests if it is present
# DEFT - default run target
# if blank then use interactive mode
# tr - run target
# for running batch jobs make uses the script $(MACH).$(tr).com from the include directory
#
# internally, tt is the variable for the test to run used in the call to the batch scripts
#

ifndef tr
   ifdef DEFT
      tr := $(DEFT)
   endif
endif
#
ifdef EIGENDIR
   IEIGEN := -I$(EIGENDIR)
else
   IEIGEN :=
endif

ifeq ($(MACH), apple)
    LIBINT_TOPDIR=/Users/yuchen/Applications/libint
endif
$(info CURDIR is $(CURDIR))
$(info pe is $(pe))
$(info MACH is $(MACH))
$(info machcomp is $(machcomp))
$(info bindir is $(bindir))
$(info objdir is $(objdir))
$(info srcdir is $(srcdir))
$(info output directory (od) is $(od))
$(info LIBINT_TOPDIR is $(LIBINT_TOPDIR))
$(info Eigen directory is $(EIGENDIR) using $(IEIGEN))
$(info git sha-1 is $(GITSHA1))
$(info git date $(GITDATE))
#
#include include/$(machcomp).mk
#
$(info RUNBATCH is $(RUNBATCH))
$(info run script is $(tr))
#
ROOTDIR := $(notdir $(pe))

# Variables used by LPP for correct preprocessing of Fortran routines
export DPRE CRAY DOSF IRIXP UAIX UAIXMC GEN LAPACKR LAPACKC MPI MPIF90 GITSHA1

VPATH = $(binrel) src $(objrel) libint_interface

SEGS_1 = Orient.o VcpPol.o ExpOrb.o

SEGS_2 = GaussianCnv.o GamessCnv.o MesaCnv.o SymGen.o GenGrid.o AngGCt.o RotOrb.o Den.o StPot.o AsyPol.o Fege.o \
   GetRMax.o ScatStabA.o SymNormMode.o SymProd.o MatEle.o DipoleOp.o GetDataRecordDef.o \
  TotalCrossSection.o MoldenCnv.o MFDCS.o EDCS.o VcpPN.o \
  VcpBN.o SaveData.o DPot.o Resonance.o ViewOrb.o \
  Vposfit.o Rotate.o OrientN.o CalcInt.o DumpMesa.o RotOrientAsym.o SchmidtOrth.o  \
   MFTimeDelay.o LFTimeDelay.o 

OVERSET_OBJS = InterpolationOverset.o G0Overset.o ExchangeOverset.o ScatStabOverset.o StPot_overset.o PatchAlgo.o

INTERFACEOBJS = 
SEGIDY = VibAve.o VibAveN.o VibAveNI.o 

SEGSO = GetPGroup.o 

SUBS1P = Sub1p.o
SUBSF90_1 = Sub.o ModuleSubs.o SlatecSpline.o
SUBSF90_2 = CrossSection.o Asymgen.o
SUBIDY = CnvIdy.o RotateIdy.o
SUBSSetUp = SubInA.o SubInB.o
SUBSLPP = SubLPP.o CoulCC.o
ifeq ($(wildcard $(objdir)/CurrSHA1.dat),)
   $(info No SHA1 found)
   $(shell echo $(GITSHA1)  > $(objdir)/CurrSHA1.dat)
   FORCELPP =  SHA1.o
else ifneq ($(shell cat $(objdir)/CurrSHA1.dat),$(GITSHA1))
   $(info New SHA1 found)
   FORCELPP =  SHA1.o
   $(shell echo $(GITSHA1) > $(objdir)/CurrSHA1.dat)
else
   $(info SHA1 not changed)
   NOFORCELPP =  SHA1.o
endif
SUBS = $(SUBSF90_1) $(SUBSF90_2) $(SUBSLPP) $(SUBSSetUp)

# LIBINT interface directory and variables
LIBINT_LIBS := -L$(LIBINT_TOPDIR)/lib -lint2 $(CPPLIBFLAGS)
LIBINT_SRCDIR  := $(pe)/libint_interface
LIBINT_OBJS := IntegrationTools.o ElectrostaticPointIntegrals.o
LIBINT_H := $(LIBINT_OBJS:.o=.h)

CPPFLAGS := -I$(LIBINT_TOPDIR)/include -I$(LIBINT_TOPDIR)/include/libint2 $(IEIGEN)
CPPFLAGS := $(CPPFLAGS) -I$(LIBINT_SRCDIR)

MODS1P = Vectors.o Modules.o SetUp.o
MODS = LMPI.o $(MODS1P)
PROGS = ePolyScat.exe Testmmp.exe
PROGS1P = dat2igor.exe dat2igor1.exe MakeGeom4igor.exe MakeManual.exe MakeManualDropPages.exe \
          CompDiff.exe TestCul.exe TestDip.exe TestDipC.exe Testmm.exe CompDiffMax.exe \
	   			ViewGeom2igor.exe FindRot.exe \
          Vib1DModel.exe CnvMathNa.exe CnvMath.exe CnvMatLabNa.exe BendOrientNa.exe BendOrient.exe \
          BendONaEnsight.exe ViewOrb2ensight.exe dat2ensight.exe MakeOverlap.exe GenFLMMP.exe Cube2igor.exe \
          Gen2DCuts.exe SumMFPAD.exe CnvMatLab.exe CnvMatLabSC.exe AddDatMatchLab.exe BendOrientSC.exe \
          PartialWave.exe BendOrientDA.exe GetPoints1D.exe GetPoints2D.exe VibAveOrth.exe CnvLinFull.exe NRFPAD.exe \
          MoldenMerge.exe dat1Dto2D.exe dat1DtoXYZ.exe datb2dat.exe Cnv2PolyDCS.exe grepdat.exe grepf90.exe \
          MoldenRotOrb.exe TestWvFun.exe \
          MFDist1Dto2D.exe datxyzto2D.exe MFDCS2ensight.exe RatioDat.exe CnvEtoP.exe CnvoEtoP.exe HARF.exe \
          MFPADFLNNu.exe SetTOrientAsym.exe AddDatMatLab.exe dat2DSwapXY.exe AddDatCoef.exe AddDatFLN.exe CnvLin1D.exe \
	   			CheckUnitT.exe datxyz2igor1.exe TestDumpIdyAll.exe
# C++/CUDA replacement for G0Overset
CPP_OBJS  := $(objdir)/G0Overset.o
CUDA_OBJS := $(objdir)/ApplyG0.o

# Combine all objects for overset
OVERSET_OBJS := $(FORT_OBJS) $(CPP_OBJS) $(CUDA_OBJS)

$(PROGS): %.exe: %.o $(MODS) $(SUBS) $(SEGS_1) $(SEGS_2) $(SEGSO) $(SEGIDY) $(SUBIDY) $(OVERSET_OBJS) ScatStabB.o ScatStabC.o SHA1.o $(LIBINT_OBJS)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -o $(bindir)/$(@) $(LIBINT_OBJS) $(LIBINT_LIBS) $(*).o $(MODS) $(SUBS) $(SEGS_1) $(SEGS_2) \
	$(SEGIDY) $(SUBIDY) $(SEGSO) $(OVERSET_OBJS) ScatStabB.o \
	ScatStabC.o SHA1.o $(LFLAGS) $(LIBINT_LIBS)

.PHONY : all
all : $(PROGS) $(PROGS1P) manual mpiwrapper permissions

.PHONY : clean
clean:
	-rm -f $(objdir)/*.o $(objdir)/*.mod $(bindir)/*.exe manual/*.html \
	$(objdir)/CurrSHA1.dat $(bindir)/ePolyScat;

.PHONY : permissions
permissions:
	chmod a+rx bin $(bindir) include manual src tests Makefile
	chmod a+r README*
	chmod a+rx $(bindir)/*.exe $(bindir)/ePolyScat
	chmod a+r bin/README obj/README
	chmod a+rx src/* tests/* include/* manual/*
	chmod a+r extras/CompDiff.sh

$(PROGS1P): %.exe: %.f90 $(SUBS1P) CoulCC.o $(MODS1P)
	cd $(objdir);\
	$(FORT) $(FFLAGS) -o $(bindir)/$(@) $(srcdir)/$(*).f90 $(MODS1P) CoulCC.o $(SUBS1P) $(LFLAGS)

$(SEGS_1): %.o: %.f90 $(MODS) $(SUBS)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

$(SEGS_2): %.o: %.f90 $(MODS) $(SUBS) $(SEGS_1)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

ePolyScat.o: ePolyScat.f90 $(MODS) $(SUBS) $(SEGS_1) $(SEGS_2) $(SEGSO) $(SEGIDY) $(OVERSET_OBJS) ScatStabB.o ScatStabC.o SHA1.o $(LIBINT_OBJS)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/ePolyScat.f90

Testmmp.o: Testmmp.f90 $(MODS) $(SUBS) $(SEGS_1) $(SEGS_2) $(SEGSO) $(SEGIDY) $(OVERSET_OBJS) ScatStabB.o ScatStabC.o SHA1.o $(LIBINT_OBJS)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/Testmmp.f90

$(LIBINT_OBJS): %.o: %.cc $(LIBINT_H)
	cd $(objdir);\
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $(objdir)/$(*).o $(LIBINT_SRCDIR)/$(*).cc

$(SEGIDY): %.o: %.f90 $(MODS) $(SUBS) $(SUBIDY)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

$(SEGSO): %.o: %.f90 $(MODS) $(SUBS)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPIO) -c $(srcdir)/$(*).f90


ScatStabB.o: %.o: %.f90 $(MODS) $(SUBS) ScatStabA.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPIO) -c $(srcdir)/$(*).f90

ScatStabC.o: %.o: %.f90 $(MODS) $(SUBS) ScatStabB.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

InterpolationOverset.o: %.o: %.f90 $(MODS) $(SUBS) $(SEGS_1)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

StPot_overset.o: %.o: %.f90 $(MODS) $(SUBS) $(SEGS_1) InterpolationOverset.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

G0Overset.o: %.o: %.f90 $(MODS) $(SUBS) $(SEGS_1) InterpolationOverset.o PatchAlgo.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

ExchangeOverset.o: %.o: %.f90 $(MODS) $(SUBS) $(SEGS_1) InterpolationOverset.o G0Overset.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

ScatStabOverset.o: %.o: %.f90 $(MODS) $(SUBS) $(SEGS_1) InterpolationOverset.o G0Overset.o ExchangeOverset.o StPot_overset.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -fallow-argument-mismatch -c $(srcdir)/$(*).f90

PatchAlgo.o: %.o: %.f90 $(MODS) $(SUBS) $(SEGS_1) InterpolationOverset.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -fallow-argument-mismatch -c $(srcdir)/$(*).f90

$(SUBSLPP) $(SUBS1P): %.o: %.for $(MODS)
	cd $(objdir);\
	$(bindir)/LPP.exe <$(srcdir)/$(*).for >$(TMPDIR)/$(*).f90;\
	$(FORTMPI) $(FFLAGSMPI) -c $(TMPDIR)/$(*).f90

$(FORCELPP) : %.o: %.for .FORCE LPP.exe
	cd $(objdir);\
	$(bindir)/LPP.exe <$(srcdir)/$(*).for >$(TMPDIR)/$(*).f90;\
	$(FORTMPI) $(FFLAGSMPI) -c $(TMPDIR)/$(*).f90

.FORCE:

$(NOFORCELPP) : %.o: %.for LPP.exe
	cd $(objdir);\
	$(bindir)/LPP.exe <$(srcdir)/$(*).for >$(TMPDIR)/$(*).f90;\
	$(FORTMPI) $(FFLAGSMPI) -c $(TMPDIR)/$(*).f90

$(SUBIDY): %.o: %.f90 $(MODS) $(SUBSF90_1) $(SUBSF90_2)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90


$(SUBSF90_1): %.o: %.f90 $(MODS)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

$(SUBSF90_2): %.o: %.f90 $(MODS) $(SUBSF90_1)
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

$(SUBSSetUp): %.o: %.f90 SetUp.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/$(*).f90

$(MODSLI):
	echo "hey dummy"

Modules.o: Modules.f90 Vectors.o SetUp.o LMPI.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/Modules.f90

Vectors.o: Vectors.m4 SetUp.o
	cd $(objdir); \
	m4 $(M4FLAGS) $(srcdir)/Vectors.m4 >$(TMPDIR)/Vectors.f90;\
	 $(FORTMPI) $(FFLAGSMPI) -c $(TMPDIR)/Vectors.f90

LMPI.o: LMPI.f90 SetUp.o
	cd $(objdir);\
	$(FORTMPI) $(FFLAGSMPI) -c $(srcdir)/LMPI.f90

SetUp.o: SetUp.for LPP.exe
	cd $(objdir);\
	$(bindir)/LPP.exe <$(srcdir)/SetUp.for >$(TMPDIR)/SetUp.f90;\
	$(FORTMPI) $(FFLAGSMPI) -c $(TMPDIR)/SetUp.f90

LPP.exe: LPP.f90
	-cd $(TMPDIR);\
	$(FORT) $(FFLAGS) -o $(bindir)/LPP.exe $(srcdir)/LPP.f90 $(LFLAGS) ;\
	rm -f *.o *.mod

testcl : 
	$(MAKE) TestCul.exe
	tests/testcl.job >&tests/testcl.out$$$$

.PHONY : mpiwrapper
mpiwrapper:
	include/$(machcomp).mpiwrapper.gen $(binrel) 

.PHONY: $(TESTS) $(TESTSX)
$(TESTS) $(TESTSX) :

ifndef tr
ifdef o
ifeq ($(o), stdout)
	$(MAKE) ePolyScat.exe;\
	export pd=$(pe)/tests/$@.d; \
	$(bindir)/ePolyScat $(pe)/tests/$@.inp
else
	$(MAKE) ePolyScat.exe;\
	export pd=$(pe)/tests/$@.d; \
	$(bindir)/ePolyScat $(pe)/tests/$@.inp >&$(o)
endif
else
	$(MAKE) ePolyScat.exe;\
	export pd=$(pe)/tests/$@.d; \
	$(bindir)/ePolyScat $(pe)/tests/$@.inp >&$(od)/$@.out
endif
else
ifeq ($(findstring batch,$(tr)), batch)
ifeq ($(RUNBATCH), bsub)
	$(MAKE) ePolyScat.exe; \
	export tt=$@ ; export e=$(pe) ; export bd=$(bindir) ; export ode=$(od) ; \
	$(pe)/include/$(MACH).$(tr).com
else
ifeq ($(RUNBATCH), sbatch)
	$(MAKE) ePolyScat.exe; \
	$(RUNBATCH) --export=tt=$@,e=$(pe),bd=$(bindir),ode=$(od),ALL $(pe)/include/$(MACH).$(tr).com
else
	$(MAKE) ePolyScat.exe; \
	$(RUNBATCH) -v tt=$@,e=$(pe),bd=$(bindir),ode=$(od) $(pe)/include/$(MACH).$(tr).com
endif
endif
else
	# bad target
endif
endif

.PHONY : testall
testall:
ifndef tr
	$(MAKE) ePolyScat.exe;\
	for testj in $(TESTS); do \
	$(bindir)/ePolyScat $(pe)/tests/$$testj.inp >& $(od)/$$testj.out ; done;
else
ifeq ($(findstring batch,$(tr)), batch)
ifeq ($(RUNBATCH), bsub)
	$(MAKE) ePolyScat.exe; \
	for testj in $(TESTS); do \
	export tt=$$testj ; export e=$(pe) ; export bd=$(bindir) ; export ode=$(od) ; \
	$(pe)/include/$(MACH).$(tr).com ; done;
else
ifeq ($(RUNBATCH), sbatch)
	$(MAKE) ePolyScat.exe; \
	for testj in $(TESTS); do \
	$(RUNBATCH) --export=tt=$$testj,e=$(pe),bd=$(bindir),ode=$(od),ALL $(pe)/include/$(MACH).$(tr).com ; done ;
else
	$(MAKE) ePolyScat.exe;\
	for testj in $(TESTS); do \
	$(RUNBATCH) -v tt=$$testj,e=$(pe),bd=$(bindir),ode=$(od) $(pe)/include/$(MACH).$(tr).com ; done;
endif
endif
else
	# bad target
endif
endif

.PHONY : testallckohn
testallckohn:
ifndef tr
	$(MAKE) ePolyScat.exe;\
	for testj in $(TESTSCKOHN); do \
	$(bindir)/ePolyScat $(pe)/tests/$$testj.inp >& $(od)/$$testj.out ; done;
else
ifeq ($(findstring batch,$(tr)), batch)
ifeq ($(RUNBATCH), bsub)
	$(MAKE) ePolyScat.exe; \
	for testj in $(TESTSCKOHN); do \
	export tt=$$testj ; export e=$(pe) ; export bd=$(bindir) ; export ode=$(od) ; \
	$(pe)/include/$(MACH).$(tr).com ; done;
else
ifeq ($(RUNBATCH), sbatch)
	$(MAKE) ePolyScat.exe; \
	for testj in $(TESTSCKOHN); do \
	$(RUNBATCH) --export=tt=$$testj,e=$(pe),bd=$(bindir),ode=$(od),ALL $(pe)/include/$(MACH).$(tr).com ; done ;
else
	$(MAKE) ePolyScat.exe;\
	for testj in $(TESTSCKOHN); do \
	$(RUNBATCH) -v tt=$$testj,e=$(pe),bd=$(bindir),ode=$(od) $(pe)/include/$(MACH).$(tr).com ; done;
endif
endif
else
	# bad target
endif
endif

.PHONY : manual
manual:
	$(MAKE) MakeManual.exe;\
	. $(pe)/manual/MakeManual.com

.PHONY : manualdroppages
manualdroppages:
	$(MAKE) MakeManualDropPages.exe;\
	. $(pe)/manualdroppages/MakeManualDropPages.com



LIBINT_TOPDIR=/Users/yuchen/Applications/libint
TOPDIR=$(LIBINT_TOPDIR)
$(info TOPDIR is $(TOPDIR))
ifndef SRCDIR
  # SRCDIR=$(shell pwd)
  # SRCDIR=$(TOPDIR)/build/
endif
# -include $(TOPDIR)/../libint/build_edison/src/bin/MakeVars
# -include $(TOPDIR)/../libint/build_cori/src/bin/MakeVars

-include $(TOPDIR)/../libint/$(LIBINT_BUILD)/MakeVars
# include headers the object include directory
# CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(TOPDIR)/../include/ $(EIGEN3)/eigen3 -I./ # -I$(SRCDIR)/src/lib/libint -DSRCDATADIR=\"$(SRCDIR)/lib/basis\"
# CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(TOPDIR)/../include/ $(IEIGEN) -I./ # -I$(SRCDIR)/src/lib/libint -DSRCDATADIR=\"$(SRCDIR)/lib/basis\"
# CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I./ $(IEIGEN) -std=c++11
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(TOPDIR)/../include/ $(IEIGEN) -I./
CPPFLAGS += -std=c++11
CPPFLAGS += -Wno-deprecated-declarations 
# CPPFLAGS += -I./eigen3 # -I$(HOME)/Projects/libint2lib_cori/include -I$(HOME)/Projects/libint2lib_cori/include/libint2
# CPPFLAGS += -I/usr/include/c++/4.7

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -lint2
vpath %.a $(TOPDIR)/lib:$(TOPDIR)/lib/.libs

OBJSUF = o
DEPSUF = d
CXXDEPENDSUF = none
CXXDEPENDFLAGS = -M


TEST1 = ElectrostaticPointIntegrals
CXXTEST1SRC = $(TEST1).cc
CXXTEST1OBJ = $(CXXTEST1SRC:%.cc=%.$(OBJSUF))
CXXTEST1DEP = $(CXXTEST1SRC:%.cc=%.$(DEPSUF))

TEST2 = IntegrationTools
CXXTEST2SRC = $(TEST2).cc
CXXTEST2OBJ = $(CXXTEST2SRC:%.cc=%.$(OBJSUF))
CXXTEST2DEP = $(CXXTEST2SRC:%.cc=%.$(DEPSUF))

interface:: $(CXXTEST1OBJ) $(CXXTEST2OBJ)
	mv -v $(CXXTEST1OBJ) $(CXXTEST2OBJ) $(objdir)

check1:: $(TEST1)

check2::  $(TEST2)

$(TEST1): $(CXXTEST1OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)

$(TEST2): $(CXXTEST2OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)

# Source files for timer and tester are to be compiled using CXXGEN
$(TEST1) $(TEST2): CXX=$(CXXGEN)
$(TEST1) $(TEST2): CXXFLAGS=$(CXXGENFLAGS)
$(TEST1) $(TEST2): LD=$(CXXGEN)

clean::
	-rm -rf $(TEST1) $(TEST2) *.o *.d

distclean:: realclean
	-rm -rf $(TOPDIR)/include/libint2/boost

realclean:: clean

targetclean:: clean

$(TOPDIR)/include/libint2/boost/preprocessor.hpp: $(SRCDIR)/external/boost.tar.gz
	gunzip -c $(SRCDIR)/external/boost.tar.gz | tar -xf - -C $(TOPDIR)/include/libint2

depend:: $(CXXTEST1DEP) $(CXXTEST2DEP)

ifneq ($(CXXDEPENDSUF),none)
%.d:: %.cc $(TOPDIR)/include/libint2/boost/preprocessor.hpp
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< > /dev/null
	sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' < $(*F).$(CXXDEPENDSUF) > $(@F)
	/bin/rm -f $(*F).$(CXXDEPENDSUF)
else
%.d:: %.cc $(TOPDIR)/include/libint2/boost/preprocessor.hpp
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< | sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' > $(@F)
endif

-include $(CXXTEST1DEP)
-include $(CXXTEST2DEP)

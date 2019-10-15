SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

VPATH += src

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
# CPLEXDIR      = /home/rottner/Programmes/Cplex12.8/cplex
# CONCERTDIR    = /home/rottner/Programmes/Cplex12.8/concert

#------------------------------------------------------------

include cplex_dir.mk

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl 
CLNFLAGS  = -lcplex -lm -lpthread -ldl 


CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include
EXDATA        = $(EXDIR)/data
EXSRCCPP      = $(EXDIR)/src/cpp

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)  -Wno-ignored-attributes
JCFLAGS = $(JOPT)


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------



# ------------------------------------------------------------
#
# Obj
#

OBJ_LIST := InstanceUCP.o \
	    ModeleUCP.o \
	    ModeleFlot.o \
	    Branching.o \
	    Node.o \
	    Mob.o \
	    Fixing.o \
        Process.o \
	    main.o \

OBJ:= $(addprefix build/, $(OBJ_LIST))



# ------------------------------------------------------------
# link

BIN      := mf

$(BIN): $(OBJ)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) $^ -o $@  $(CCLNFLAGS)
	#$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $@ $^ $(CCLNFLAGS)


# ------------------------------------------------------------
# compil

build/%.o: src/%.cc
	$(CCC) -c $(CCFLAGS) $< -o $@


# ------------------------------------------------------------
.PHONY: clean mrproper doc

clean:
	@$(remove) $(OBJ)

mrproper: clean
	@$(remove) $(BIN)



# List of dependencies
build/InstanceUCP.o   : InstanceUCP.h
build/ModeleFlot.o : InstanceUCP.h ModeleFlot.h
build/ModeleUCP.o : InstanceUCP.h ModeleUCP.h
build/Branching.o : Node.h InstanceUCP.h
build/Node.o : Node.h InstanceUCP.h
build/Mob.o : Mob.h Node.h
build/Fixing.o : Fixing.h Node.h
build/Process.o : Process.h InstanceUCP.h ModeleUCP.h
build/main.o      : InstanceUCP.h ModeleUCP.h Node.h Fixing.h Process.h Mob.h

# Local Variables:
# mode: makefile
# End:

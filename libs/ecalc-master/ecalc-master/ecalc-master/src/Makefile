# Change these for where the bioplib includes and libraries are stored
INCFLAGS = -I$(HOME)/include
LIBFLAGS = -L$(HOME)/lib

# Change these for your compiler and linker
CC      = cc $(INCFLAGS) -c -O3 -ansi -pedantic -Wall
LINK    = cc $(LIBFLAGS) -o $(EXEFILE) -O3

# Change this to specify the required libraries
LIBRARIES = -lbiop -lgen -lm -lxml2

# Change this for your delete command
DEL = \rm -f

# Change this for your prototype generator
PROTO = $(HOME)/bin/ansi -p

OFILES  = main.o ReadParams.o PrintParams.o ReadRTop.o BuildTop.o \
GetParam.o energy.o ReadStructure.o shake.o

PFILES = BuildTop.p GetParam.p energy.p ReadStructure.p shake.p

EXEFILE = ecalc

$(EXEFILE) : $(OFILES)
	$(LINK) $(OFILES) $(LIBRARIES)



main.o : main.c ecalc.h params.h
	$(CC) main.c

ReadParams.o : ReadParams.c ecalc.h params.h
	$(CC) ReadParams.c

PrintParams.o : PrintParams.c params.h
	$(CC) PrintParams.c

ReadRTop.o : ReadRTop.c ecalc.h rtop.h
	$(CC) ReadRTop.c

BuildTop.o : BuildTop.c rtop.h params.h topology.h ecalc.h
	$(CC) BuildTop.c

GetParam.o : GetParam.c params.h ecalc.h
	$(CC) GetParam.c

energy.o : energy.c topology.h ecalc.h
	$(CC) energy.c

ReadStructure.o : ReadStructure.c topology.h ecalc.h
	$(CC) ReadStructure.c


shake.o : shake.c topology.h ecalc.h
	$(CC) shake.c

protos : $(PFILES)


ecalc.h : topology.h params.h protos.h
	touch ecalc.h

.c.p : 
	$(PROTO) $< $@

clean :
	$(DEL) $(OFILES)

distclean : clean
	$(DEL) $(EXEFILE)

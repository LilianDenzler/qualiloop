# You are advised to change this to provide warning messages for your
# compiler - the sample line is for GCC. If you are compiling for
CC    = cc -O3
#CC    = cc -O3 -ansi -Wall -pedantic 
# Mac OS/X, you may need the full path to gcc
#CC    = /Developer/SDKs/MacOSX10.4u.sdk/usr/bin/gcc -O3

# Mac OS/X, uncomment the following
#MAC = -mmacosx-version-min=10.4 -L/Developer/SDKs/MacOSX10.4u.sdk/usr/lib

# To allow use of the GNU Readline library, uncomment the following two
# lines - You may need to install the GNU readline development libraries
# first!
#READLINE    = -DREADLINE_SUPPORT
#READLINELIB = -lreadline -lcurses

# Uncomment if you want to use the rotate and refit code for avoiding 
# local minimum problem. If compiled with gcc -O3 this problem seems to
# go away in any case and only affects fitting of identical structures
# A fix will soon be provided in fit.c instead
#ROTATEREFIT = -DROTATE_REFIT

# To allow decompression of gzipped PDB files on the fly, change the 
# COPT line to
COPT  = -DGUNZIP_SUPPORT


# Uncomment these and fix the paths if you have the XMAS library and
# want XMAS support
# XMAS  = -DUSE_XMAS  -I/home/amartin/inpharmatica/cvs/software/include
# XMASLIB = -L/home/amartin/inpharmatica/cvs/software/lib -lxmas -lacrm

LINK1 = 
LINK2 =
ANSI  = ansi -p 

OFILES = main.o todo.o fitting.o NWAlign.o
PFILES = main.p todo.p fitting.p NWAlign.p
IFILES = ProFit.h protos.h
LFILES = bioplib/WritePDB.o bioplib/parse.o \
bioplib/PDB2Seq.o bioplib/WindIO.o bioplib/help.o \
bioplib/fsscanf.o bioplib/ApMatPDB.o bioplib/align.o \
bioplib/fit.o bioplib/throne.o bioplib/angle.o bioplib/array2.o \
bioplib/ParseRes.o bioplib/OpenFile.o bioplib/padterm.o \
bioplib/ReadPIR.o bioplib/openorpipe.o bioplib/GetWord.o \
bioplib/SelAtPDB.o bioplib/IndexPDB.o bioplib/CopyPDB.o \
bioplib/FindZonePDB.o bioplib/countchar.o bioplib/KillLeadSpaces.o \
bioplib/MatMult3_33.o bioplib/StringToUpper.o bioplib/upstrncmp.o \
bioplib/LegalAtomSpec.o bioplib/GetPDBChainLabels.o bioplib/DupePDB.o \
bioplib/TranslatePDB.o bioplib/AtomNameMatch.o bioplib/chindex.o \
bioplib/FindNextResidue.o bioplib/StoreString.o bioplib/FreeStringList.o \
bioplib/MatMult33_33.o bioplib/ReadPDB.o bioplib/aalist.o \
bioplib/CreateRotMat.o bioplib/BuildConect.o bioplib/hash.o \
bioplib/stringutil.o bioplib/prime.o bioplib/FindResidue.o

profit : $(OFILES) $(LFILES)
	$(CC) -o profit $(MAC) $(OFILES) $(LFILES) $(XMAS) -lm $(LINK2) $(XMASLIB) $(READLINELIB)

.c.o : $(IFILES)
	$(CC) $(COPT) $(MAC) $(READLINE) $(ROTATEREFIT) -o $@ -c $<

clean :
	/bin/rm $(OFILES) $(LFILES)

protos : $(PFILES)

.c.p   :
	$(ANSI) $< $@

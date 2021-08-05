CC      = gcc
BINDIR  = /serv/www/html_lilian/libs//bin
DATADIR = /serv/www/html_lilian/libs//data
CFLAGS  = -O3 -ansi -Wall -pedantic -Wno-stringop-truncation -I./bioplib -L./bioplib
LFLAGS  = -lbiop -lgen -lm -lxml2
TARGETS = scorecons pdbgetzone pdbrenum pdbavbr naccess2bval pdbdummystrip pdbconect pdbsplitchains pdbmakepatch pdborder pdbtranslate pdbhadd pdb2ms pdbhstrip rangecontacts pdbatomcount pdbfindresrange pdbtorsions setpdbnumbering pdb2pir sixft pdbflip pdb2pdbml pdbhphob pdbatomselect pdbaddhet distmat pdbsymm chaincontacts pdbpatchbval pdbcalcrms pdbml2pdb pdblistss pdbgetresidues pdbcter pdbheader pdbfit pdbatoms pdbchain pdbsphere pdbfindnearres pdbcheckforres pdbpatchnumbering pdbcentralres pdbsecstr pdbsumbval pdbgetchain pdbline pdbhbond pdbrotate pdbselect pdbcount checkpdb pdborigin pdbsolv pdb2xyz pdbhetstrip 

all : bioplib $(TARGETS)
	@echo " "
	@echo " --- BUILD COMPLETE --- "
	@echo " "
	@echo "To use some programs, you will need to set the environment variable"
	@echo "DATADIR to $(DATADIR)"

bioplib :
	(cd libsrc/bioplib/src; make)
	mkdir -p bioplib
	cp libsrc/bioplib/src/*.[ah] bioplib
	if [ ! -d ../data ]; then mkdir ../data; fi
	cp -R libsrc/bioplib/data/* ../data


install : 
	mkdir -p $(BINDIR)
	cp $(TARGETS) $(BINDIR)
	mkdir -p $(DATADIR)
	./installdata.sh $(DATADIR)
	@echo " "
	@echo " --- INSTALL COMPLETE --- "
	@echo " "
	@echo "To use some programs, you will need to set the environment variable"
	@echo "DATADIR to $(DATADIR)"


distclean : clean
	\rm Makefile

clean : 
	\rm -rf bioplib
	(cd libsrc/bioplib/src; make clean)
	\rm -f $(TARGETS)


links : 
	(cd $(BINDIR); \rm -f addhet      ; ln -s pdbaddhet          addhet       )
	(cd $(BINDIR); \rm -f as2bval     ; ln -s naccess2bval       as2bval      )
	(cd $(BINDIR); \rm -f atomcount   ; ln -s pdbatomcount       atomcount    )
	(cd $(BINDIR); \rm -f avbr        ; ln -s pdbavbr            avbr         )
	(cd $(BINDIR); \rm -f centralres  ; ln -s pdbcentralres      centralres   )
	(cd $(BINDIR); \rm -f chainpdb    ; ln -s pdbchain           chainpdb     )
	(cd $(BINDIR); \rm -f checkforres ; ln -s pdbcheckforres     checkforres  )
	(cd $(BINDIR); \rm -f countpdb    ; ln -s pdbcount           countpdb     )
	(cd $(BINDIR); \rm -f dummystrip  ; ln -s pdbdummystrip      dummystrip   )
	(cd $(BINDIR); \rm -f findresrange; ln -s pdbfindresrange    findresrange )
	(cd $(BINDIR); \rm -f flip        ; ln -s pdbflip            flip         )
	(cd $(BINDIR); \rm -f getchain    ; ln -s pdbgetchain        getchain     )
	(cd $(BINDIR); \rm -f getpdb      ; ln -s pdbgetzone         getpdb       )
	(cd $(BINDIR); \rm -f getresidues ; ln -s pdbgetresidues     getresidues  )
	(cd $(BINDIR); \rm -f hetstrip    ; ln -s pdbhetstrip        hetstrip     )
	(cd $(BINDIR); \rm -f hstrip      ; ln -s pdbhstrip          hstrip       )
	(cd $(BINDIR); \rm -f makepatch   ; ln -s pdbmakepatch       makepatch    )
	(cd $(BINDIR); \rm -f nullstrip   ; ln -s pdbdummystrip      nullstrip    )
	(cd $(BINDIR); \rm -f numpdb      ; ln -s setpdbnumbering    numpdb       )
	(cd $(BINDIR); \rm -f patchbval   ; ln -s pdbpatchbval       patchbval    )
	(cd $(BINDIR); \rm -f patchpdbnum ; ln -s pdbpatchnumbering  patchpdbnum  )
	(cd $(BINDIR); \rm -f renumpdb    ; ln -s pdbrenum           renumpdb     )
	(cd $(BINDIR); \rm -f rmspdb      ; ln -s pdbcalcrms         rmspdb       )
	(cd $(BINDIR); \rm -f rotate      ; ln -s pdbrotate          rotate       )
	(cd $(BINDIR); \rm -f solv        ; ln -s pdbsolv            solv         )
	(cd $(BINDIR); \rm -f splitchains ; ln -s pdbsplitchains     splitchains  )
	(cd $(BINDIR); \rm -f sumbval     ; ln -s pdbsumbval         sumbval      )
	(cd $(BINDIR); \rm -f torsions    ; ln -s pdbtorsions        torsions     )
	(cd $(BINDIR); \rm -f transpdb    ; ln -s pdbtranslate       transpdb     )
	(cd $(BINDIR); \rm -f pdbatomsel  ; ln -s pdbatomselect      pdbatomsel   )
	(cd $(BINDIR); \rm -f atomsel     ; ln -s pdbatomselect      atomsel      )


scorecons : scorecons.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbgetzone : pdbgetzone.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbrenum : pdbrenum.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbavbr : pdbavbr.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

naccess2bval : naccess2bval.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbdummystrip : pdbdummystrip.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbconect : pdbconect.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbsplitchains : pdbsplitchains.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbmakepatch : pdbmakepatch.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdborder : pdborder.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbtranslate : pdbtranslate.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbhadd : pdbhadd.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdb2ms : pdb2ms.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbhstrip : pdbhstrip.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

rangecontacts : rangecontacts.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbatomcount : pdbatomcount.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbfindresrange : pdbfindresrange.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbtorsions : pdbtorsions.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

setpdbnumbering : setpdbnumbering.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdb2pir : pdb2pir.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

sixft : sixft.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbflip : pdbflip.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdb2pdbml : pdb2pdbml.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbhphob : pdbhphob.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbatomselect : pdbatomselect.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbaddhet : pdbaddhet.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

distmat : distmat.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbsymm : pdbsymm.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

chaincontacts : chaincontacts.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbpatchbval : pdbpatchbval.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbcalcrms : pdbcalcrms.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbml2pdb : pdbml2pdb.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdblistss : pdblistss.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbgetresidues : pdbgetresidues.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbcter : pdbcter.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbheader : pdbheader.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbfit : pdbfit.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbatoms : pdbatoms.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbchain : pdbchain.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbsphere : pdbsphere.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbfindnearres : pdbfindnearres.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbcheckforres : pdbcheckforres.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbpatchnumbering : pdbpatchnumbering.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbcentralres : pdbcentralres.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbsecstr : pdbsecstr.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbsumbval : pdbsumbval.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbgetchain : pdbgetchain.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbline : pdbline.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbhbond : pdbhbond.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbrotate : pdbrotate.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbselect : pdbselect.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbcount : pdbcount.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

checkpdb : checkpdb.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdborigin : pdborigin.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbsolv : pdbsolv.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdb2xyz : pdb2xyz.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

pdbhetstrip : pdbhetstrip.c
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

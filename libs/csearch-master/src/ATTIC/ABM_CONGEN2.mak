ROOT      = $(OML_ABM)
CINCPATH  = -I$(ROOT)/INC -I$(ROOT)/PROTO
FINCPATH  = -I$(ROOT)/INC
OBJDIR    = $(ROOT)/OBJ
LIBDIR    = $(ROOT)/LIB
CFLAGS    = -prototypes -acpp -G 4 -D__SILICONGRAPHICS__ -D__IRIX__ -D__IRIX_3_3__ -D__FORT_UNDERSCORE__ -D__UNIX__ -DSYSV -DXOPEN_CATALOG -DOUX_GL -cckr
FFLAGS    = -Nc30
 
 
.f.a:
	f77 -c $(FFLAGS) $(FINCPATH) $*.f
 
.c.a:
	cc -c $(CFLAGS) $(CINCPATH) $*.c
 
 
LIB = $(LIBDIR)/ABM_CONGEN.a
 
.PRECIOUS:	$(LIB)
 
$(LIB): \
	$(LIB)(abmpad.o) \
	$(LIB)(indexf.o) \
	$(LIB)(indexi.o) \
	$(LIB)(infrls.o) \
	$(LIB)(init_debug.o) \
	$(LIB)(initfcgen.o) \
	$(LIB)(initfgener.o) \
	$(LIB)(interpa.o) \
	$(LIB)(isctrnglst.o) \
	$(LIB)(killspcs.o) \
	$(LIB)(ljust.o) \
	$(LIB)(ljustpad.o) \
	$(LIB)(lookupname.o) \
	$(LIB)(makeh.o) \
	$(LIB)(makinb.o) \
	$(LIB)(mapgrdtorng.o) \
	$(LIB)(match.o) \
	$(LIB)(matom.o) \
	$(LIB)(maxcoor.o) \
	$(LIB)(maxf2.o) \
	$(LIB)(maxi2.o) \
	$(LIB)(maxr.o) \
	$(LIB)(merfi.o) \
	$(LIB)(mincoor.o) \
	$(LIB)(minf2.o) \
	$(LIB)(minf3.o) \
	$(LIB)(mini2.o) \
	$(LIB)(mkprolrng.o) \
	$(LIB)(mygetparbond.o) \
	$(LIB)(mynextwd.o) \
	$(LIB)(nindx.o) \
	$(LIB)(onethr.o) \
	$(LIB)(pack2.o) \
	$(LIB)(pad.o) \
	$(LIB)(padspace.o) \
	$(LIB)(padterm.o) \
	$(LIB)(parse.o) \
	$(LIB)(patexh.o) \
	$(LIB)(phi.o) \
	$(LIB)(phia.o) \
	$(LIB)(prdie.o) \
	$(LIB)(print_atom.o) \
	$(LIB)(print_atom1.o) \
	$(LIB)(print_atoms.o) \
	$(LIB)(procang.o) \
	$(LIB)(processbonds.o) \
	$(LIB)(prochbnd.o) \
	$(LIB)(procimps.o) \
	$(LIB)(procnbnds.o) \
	$(LIB)(proctors.o) \
	$(LIB)(propbbones.o) \
	$(LIB)(prtcgncmnd.o) \
	$(LIB)(prtdbgvars.o) \
	$(LIB)(prtglobopt.o) \
	$(LIB)(prtgrdsz.o) \
	$(LIB)(prtrnglst.o) \
	$(LIB)(prtsidchn.o) \
	$(LIB)(pscgncmnd.o) \
	$(LIB)(rbest_f.o) \
	$(LIB)(rdemap.o) \
	$(LIB)(rdpepmaps.o) \
	$(LIB)(rdprocons.o) \
	$(LIB)(rdprtitle.o) \
	$(LIB)(readcoords.o) \
	$(LIB)(readparams.o) \
	$(LIB)(readpdb.o) \
	$(LIB)(readpgp.o) \
	$(LIB)(readpir.o) \
	$(LIB)(readrtf.o) \
	$(LIB)(regrid.o) \
	$(LIB)(rnhbsrch.o) \
	$(LIB)(rnnbsrch.o) \
	$(LIB)(rtfacceptor.o) \
	$(LIB)(rtfangle.o) \
	$(LIB)(rtfatom.o) \
	$(LIB)(rtfbond.o) \
	$(LIB)(rtfbuild.o) \
	$(LIB)(rtfdeclare.o) \
	$(LIB)(rtfdonor.o) \
	$(LIB)(rtfgroup.o) \
	$(LIB)(rtfimproper.o) \
	$(LIB)(rtfprint.o) \
	$(LIB)(rtfresidue.o) \
	$(LIB)(rtftopology.o) \
	$(LIB)(rtftorsion.o) \
	$(LIB)(rtftype.o) \
	$(LIB)(schdls.o) \
	$(LIB)(schnatm.o) \
	$(LIB)(scnemapln.o) \
	$(LIB)(selsidetor.o) \
	$(LIB)(setavdclp.o) \
	$(LIB)(setchncls.o) \
	$(LIB)(setdefs.o) \
	$(LIB)(setrcntcts.o) \
	$(LIB)(setrestitr.o) \
	$(LIB)(setsdbound.o) \
	$(LIB)(setsddofrs.o) \
	$(LIB)(setsidchn.o) \
	$(LIB)(setsidphis.o) \
	$(LIB)(setup_consp.o) \
	$(LIB)(setup_hbond.o) \
	$(LIB)(setup_imove.o) \
	$(LIB)(setup_nbond.o) \
	$(LIB)(setup_parmno.o) \
	$(LIB)(setup_parser.o) \
	$(LIB)(setup_status.o) \
	$(LIB)(setupseq.o) \
	$(LIB)(shufcombs.o) \
	$(LIB)(sidchnall.o) \
	$(LIB)(sidchnbst.o) \
	$(LIB)(sidchncomb.o) \
	$(LIB)(sidchnengy.o) \
	$(LIB)(sidchnfst.o) \
	$(LIB)(sidchnindi.o) \
	$(LIB)(sidchnitr.o) \
	$(LIB)(sidchnrms.o) \
	$(LIB)(sidechain_f.o) \
	$(LIB)(sign.o) \
	$(LIB)(skiptitle.o) \
	$(LIB)(sort_i.o) \
	$(LIB)(sort_ui.o) \
	$(LIB)(sorti_perm.o) \
	$(LIB)(srchint.o) \
	$(LIB)(srchnat1.o) \
	$(LIB)(srchnatm1.o) \
	$(LIB)(srchwd.o) \
	$(LIB)(srwdbd.o) \
	$(LIB)(statprn.o) \
	$(LIB)(status_f.o) \
	$(LIB)(stread.o) \
	$(LIB)(struppr.o) \
	$(LIB)(stuphb1.o) \
	$(LIB)(termcgen.o) \
	$(LIB)(terminate.o) \
	$(LIB)(trace.o) \
	$(LIB)(trimic.o) \
	$(LIB)(typeinres.o) \
	$(LIB)(unpack2.o) \
	$(LIB)(updtnictot.o) \
	$(LIB)(user_cgeval.o) \
	$(LIB)(vlen.o) \
	$(LIB)(write_status.o) \
	$(LIB)(writepdb.o) \
	$(LIB)(wrtcoordsf.o) \
	$(LIB)(wrtpdbrec.o) \
	$(LIB)(xbyind.o)  
	ar rcv $(LIB) *.o
	/bin/rm *.o
	@echo Completed build of ABM_CONGEN2
 
 
$(LIB)(abmpad.o):	abmpad.c
 
$(LIB)(indexf.o):	indexf.c
 
$(LIB)(indexi.o):	indexi.c
 
$(LIB)(infrls.o):	infrls.f
 
$(LIB)(init_debug.o):	init_debug.c
 
$(LIB)(initfcgen.o):	initfcgen.c
 
$(LIB)(initfgener.o):	initfgener.c
 
$(LIB)(interpa.o):	interpa.f
 
$(LIB)(isctrnglst.o):	isctrnglst.c
 
$(LIB)(killspcs.o):	killspcs.c
 
$(LIB)(ljust.o):	ljust.c
 
$(LIB)(ljustpad.o):	ljustpad.c
 
$(LIB)(lookupname.o):	lookupname.c
 
$(LIB)(makeh.o):	makeh.c
 
$(LIB)(makinb.o):	makinb.c
 
$(LIB)(mapgrdtorng.o):	mapgrdtorng.c
 
$(LIB)(match.o):	match.c
 
$(LIB)(matom.o):	matom.c
 
$(LIB)(maxcoor.o):	maxcoor.c
 
$(LIB)(maxf2.o):	maxf2.c
 
$(LIB)(maxi2.o):	maxi2.c
 
$(LIB)(maxr.o):	maxr.f
 
$(LIB)(merfi.o):	merfi.f
 
$(LIB)(mincoor.o):	mincoor.c
 
$(LIB)(minf2.o):	minf2.c
 
$(LIB)(minf3.o):	minf3.c
 
$(LIB)(mini2.o):	mini2.c
 
$(LIB)(mkprolrng.o):	mkprolrng.f
 
$(LIB)(mygetparbond.o):	mygetparbond.c
 
$(LIB)(mynextwd.o):	mynextwd.c
 
$(LIB)(nindx.o):	nindx.f
 
$(LIB)(onethr.o):	onethr.c
 
$(LIB)(pack2.o):	pack2.c
 
$(LIB)(pad.o):	pad.c
 
$(LIB)(padspace.o):	padspace.c
 
$(LIB)(padterm.o):	padterm.c
 
$(LIB)(parse.o):	parse.c
 
$(LIB)(patexh.o):	patexh.c
 
$(LIB)(phi.o):	phi.c
 
$(LIB)(phia.o):	phia.c
 
$(LIB)(prdie.o):	prdie.c
 
$(LIB)(print_atom.o):	print_atom.c
 
$(LIB)(print_atom1.o):	print_atom1.c
 
$(LIB)(print_atoms.o):	print_atoms.c
 
$(LIB)(procang.o):	procang.c
 
$(LIB)(processbonds.o):	processbonds.c
 
$(LIB)(prochbnd.o):	prochbnd.c
 
$(LIB)(procimps.o):	procimps.c
 
$(LIB)(procnbnds.o):	procnbnds.c
 
$(LIB)(proctors.o):	proctors.c
 
$(LIB)(propbbones.o):	propbbones.c
 
$(LIB)(prtcgncmnd.o):	prtcgncmnd.c
 
$(LIB)(prtdbgvars.o):	prtdbgvars.c
 
$(LIB)(prtglobopt.o):	prtglobopt.c
 
$(LIB)(prtgrdsz.o):	prtgrdsz.c
 
$(LIB)(prtrnglst.o):	prtrnglst.c
 
$(LIB)(prtsidchn.o):	prtsidchn.c
 
$(LIB)(pscgncmnd.o):	pscgncmnd.c
 
$(LIB)(rbest_f.o):	rbest_f.c
 
$(LIB)(rdemap.o):	rdemap.c
 
$(LIB)(rdpepmaps.o):	rdpepmaps.c
 
$(LIB)(rdprocons.o):	rdprocons.c
 
$(LIB)(rdprtitle.o):	rdprtitle.c
 
$(LIB)(readcoords.o):	readcoords.c
 
$(LIB)(readparams.o):	readparams.c
 
$(LIB)(readpdb.o):	readpdb.c
 
$(LIB)(readpgp.o):	readpgp.c
 
$(LIB)(readpir.o):	readpir.c
 
$(LIB)(readrtf.o):	readrtf.c
 
$(LIB)(regrid.o):	regrid.c
 
$(LIB)(rnhbsrch.o):	rnhbsrch.f
 
$(LIB)(rnnbsrch.o):	rnnbsrch.f
 
$(LIB)(rtfacceptor.o):	rtfacceptor.c
 
$(LIB)(rtfangle.o):	rtfangle.c
 
$(LIB)(rtfatom.o):	rtfatom.c
 
$(LIB)(rtfbond.o):	rtfbond.c
 
$(LIB)(rtfbuild.o):	rtfbuild.c
 
$(LIB)(rtfdeclare.o):	rtfdeclare.c
 
$(LIB)(rtfdonor.o):	rtfdonor.c
 
$(LIB)(rtfgroup.o):	rtfgroup.c
 
$(LIB)(rtfimproper.o):	rtfimproper.c
 
$(LIB)(rtfprint.o):	rtfprint.c
 
$(LIB)(rtfresidue.o):	rtfresidue.c
 
$(LIB)(rtftopology.o):	rtftopology.c
 
$(LIB)(rtftorsion.o):	rtftorsion.c
 
$(LIB)(rtftype.o):	rtftype.c
 
$(LIB)(schdls.o):	schdls.f
 
$(LIB)(schnatm.o):	schnatm.c
 
$(LIB)(scnemapln.o):	scnemapln.c
 
$(LIB)(selsidetor.o):	selsidetor.c
 
$(LIB)(setavdclp.o):	setavdclp.c
 
$(LIB)(setchncls.o):	setchncls.c
 
$(LIB)(setdefs.o):	setdefs.c
 
$(LIB)(setrcntcts.o):	setrcntcts.c
 
$(LIB)(setrestitr.o):	setrestitr.c
 
$(LIB)(setsdbound.o):	setsdbound.c
 
$(LIB)(setsddofrs.o):	setsddofrs.c
 
$(LIB)(setsidchn.o):	setsidchn.c
 
$(LIB)(setsidphis.o):	setsidphis.c
 
$(LIB)(setup_consp.o):	setup_consp.c
 
$(LIB)(setup_hbond.o):	setup_hbond.c
 
$(LIB)(setup_imove.o):	setup_imove.c
 
$(LIB)(setup_nbond.o):	setup_nbond.c
 
$(LIB)(setup_parmno.o):	setup_parmno.c
 
$(LIB)(setup_parser.o):	setup_parser.c
 
$(LIB)(setup_status.o):	setup_status.c
 
$(LIB)(setupseq.o):	setupseq.c
 
$(LIB)(shufcombs.o):	shufcombs.c
 
$(LIB)(sidchnall.o):	sidchnall.c
 
$(LIB)(sidchnbst.o):	sidchnbst.c
 
$(LIB)(sidchncomb.o):	sidchncomb.c
 
$(LIB)(sidchnengy.o):	sidchnengy.c
 
$(LIB)(sidchnfst.o):	sidchnfst.c
 
$(LIB)(sidchnindi.o):	sidchnindi.c
 
$(LIB)(sidchnitr.o):	sidchnitr.c
 
$(LIB)(sidchnrms.o):	sidchnrms.c
 
$(LIB)(sidechain_f.o):	sidechain_f.c
 
$(LIB)(sign.o):	sign.c
 
$(LIB)(skiptitle.o):	skiptitle.c
 
$(LIB)(sort_i.o):	sort_i.c
 
$(LIB)(sort_ui.o):	sort_ui.c
 
$(LIB)(sorti_perm.o):	sorti_perm.c
 
$(LIB)(srchint.o):	srchint.c
 
$(LIB)(srchnat1.o):	srchnat1.f
 
$(LIB)(srchnatm1.o):	srchnatm1.f
 
$(LIB)(srchwd.o):	srchwd.c
 
$(LIB)(srwdbd.o):	srwdbd.c
 
$(LIB)(statprn.o):	statprn.c
 
$(LIB)(status_f.o):	status_f.c
 
$(LIB)(stread.o):	stread.c
 
$(LIB)(struppr.o):	struppr.c
 
$(LIB)(stuphb1.o):	stuphb1.f
 
$(LIB)(termcgen.o):	termcgen.c
 
$(LIB)(terminate.o):	terminate.c
 
$(LIB)(trace.o):	trace.c
 
$(LIB)(trimic.o):	trimic.c
 
$(LIB)(typeinres.o):	typeinres.c
 
$(LIB)(unpack2.o):	unpack2.c
 
$(LIB)(updtnictot.o):	updtnictot.c
 
$(LIB)(user_cgeval.o):	user_cgeval.c
 
$(LIB)(vlen.o):	vlen.c
 
$(LIB)(write_status.o):	write_status.c
 
$(LIB)(writepdb.o):	writepdb.c
 
$(LIB)(wrtcoordsf.o):	wrtcoordsf.c
 
$(LIB)(wrtpdbrec.o):	wrtpdbrec.c
 
$(LIB)(xbyind.o):	xbyind.c

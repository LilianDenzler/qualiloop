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
	$(LIB)(abmerf.o) \
	$(LIB)(abmerfc.o) \
	$(LIB)(abmerfi.o) \
	$(LIB)(adafls.o) \
	$(LIB)(adatmstgrd.o) \
	$(LIB)(adatmtgrd.o) \
	$(LIB)(adccangs.o) \
	$(LIB)(addh.o) \
	$(LIB)(alloc.o) \
	$(LIB)(allstk.o) \
	$(LIB)(angle.o) \
	$(LIB)(ascale.o) \
	$(LIB)(assign_name.o) \
	$(LIB)(atomangle.o) \
	$(LIB)(backbone_f.o) \
	$(LIB)(batch_job.o) \
	$(LIB)(bin_search.o) \
	$(LIB)(binschr.o) \
	$(LIB)(bondl.o) \
	$(LIB)(c_print_atom.o) \
	$(LIB)(cartx2.o) \
	$(LIB)(cclsngdst.o) \
	$(LIB)(ccterxyz.o) \
	$(LIB)(cgen.o) \
	$(LIB)(cgrdconsp.o) \
	$(LIB)(cgrdinit.o) \
	$(LIB)(cgrdncons.o) \
	$(LIB)(cgrdnxt.o) \
	$(LIB)(cgwrtinit.o) \
	$(LIB)(cgwrtnxt.o) \
	$(LIB)(checknatc.o) \
	$(LIB)(checkrtf.o) \
	$(LIB)(chkbbone.o) \
	$(LIB)(chkchnclsr.o) \
	$(LIB)(chkclsdst.o) \
	$(LIB)(chkcntcts.o) \
	$(LIB)(chkftres.o) \
	$(LIB)(chmceil.o) \
	$(LIB)(chnclsr_f.o) \
	$(LIB)(clear_atom.o) \
	$(LIB)(clear_atoms.o) \
	$(LIB)(clear_side.o) \
	$(LIB)(clschna.o) \
	$(LIB)(cnsbblist.o) \
	$(LIB)(cnsfblist.o) \
	$(LIB)(codes.o) \
	$(LIB)(compspcgrd.o) \
	$(LIB)(constatm.o) \
	$(LIB)(constclmp.o) \
	$(LIB)(copyfarray.o) \
	$(LIB)(copyst.o) \
	$(LIB)(cpfbest.o) \
	$(LIB)(cpftbest.o) \
	$(LIB)(cpintobst.o) \
	$(LIB)(cprint.o) \
	$(LIB)(cprint1.o) \
	$(LIB)(cpyprm.o) \
	$(LIB)(csearch.o) \
	$(LIB)(debug.o) \
	$(LIB)(delatmfgrd.o) \
	$(LIB)(delatmfgrd1.o) \
	$(LIB)(delatmsfgrd.o) \
	$(LIB)(die.o) \
	$(LIB)(dispatch.o) \
	$(LIB)(dot.o) \
	$(LIB)(dststop.o) \
	$(LIB)(dump_grid.o) \
	$(LIB)(dwprangl.o) \
	$(LIB)(ecchkcont.o) \
	$(LIB)(ecntrl.o) \
	$(LIB)(ecsetlims.o) \
	$(LIB)(ecsetxyzsys.o) \
	$(LIB)(ephi.o) \
	$(LIB)(evaluate_f.o) \
	$(LIB)(explcntcts.o) \
	$(LIB)(extnlims.o) \
	$(LIB)(fill2.o) \
	$(LIB)(fill4.o) \
	$(LIB)(fill_grid.o) \
	$(LIB)(fill_grid1.o) \
	$(LIB)(fillatm_c.o) \
	$(LIB)(fillatm_ca.o) \
	$(LIB)(fillatm_cb.o) \
	$(LIB)(fillatm_h.o) \
	$(LIB)(fillatm_ht1.o) \
	$(LIB)(fillatm_n.o) \
	$(LIB)(fillatm_o.o) \
	$(LIB)(fillbbatms.o) \
	$(LIB)(fillbndang.o) \
	$(LIB)(filllog.o) \
	$(LIB)(filspc.o) \
	$(LIB)(finbbone.o) \
	$(LIB)(finchncls.o) \
	$(LIB)(finclsatm.o) \
	$(LIB)(finres.o) \
	$(LIB)(finsidchn.o) \
	$(LIB)(fixctero.o) \
	$(LIB)(fixinitgrd.o) \
	$(LIB)(fixnterh.o) \
	$(LIB)(fixpdb.o) \
	$(LIB)(forprtatm.o) \
	$(LIB)(forscanf.o) \
	$(LIB)(fortrace.o) \
	$(LIB)(free_range.o) \
	$(LIB)(fresrk.o) \
	$(LIB)(fwerf.o) \
	$(LIB)(gammln.o) \
	$(LIB)(gammp.o) \
	$(LIB)(gammq.o) \
	$(LIB)(gcd.o) \
	$(LIB)(gcf.o) \
	$(LIB)(generate.o) \
	$(LIB)(genh.o) \
	$(LIB)(genic.o) \
	$(LIB)(get_atnum.o) \
	$(LIB)(get_resnum.o) \
	$(LIB)(getclsangle.o) \
	$(LIB)(getclsbond.o) \
	$(LIB)(getdate.o) \
	$(LIB)(getgrdspc.o) \
	$(LIB)(getparam.o) \
	$(LIB)(getparbond.o) \
	$(LIB)(getparbond2.o) \
	$(LIB)(getres.o) \
	$(LIB)(getseg.o) \
	$(LIB)(getstring.o) \
	$(LIB)(gettime.o) \
	$(LIB)(getuv.o) \
	$(LIB)(gser.o) \
	$(LIB)(gtclschngom.o) \
	$(LIB)(gtprangclsa.o) \
	$(LIB)(gtprangl.o) \
	$(LIB)(gtprangl2.o) \
	$(LIB)(gtrbstatms.o) \
	$(LIB)(hadd.o) \
	$(LIB)(hbenergy.o)  
	ar rcv $(LIB) *.o
	/bin/rm *.o
	@echo Completed build of ABM_CONGEN
 
 
$(LIB)(abmerf.o):	abmerf.f
 
$(LIB)(abmerfc.o):	abmerfc.f
 
$(LIB)(abmerfi.o):	abmerfi.f
 
$(LIB)(adafls.o):	adafls.c
 
$(LIB)(adatmstgrd.o):	adatmstgrd.c
 
$(LIB)(adatmtgrd.o):	adatmtgrd.c
 
$(LIB)(adccangs.o):	adccangs.c
 
$(LIB)(addh.o):	addh.c
 
$(LIB)(alloc.o):	alloc.c
 
$(LIB)(allstk.o):	allstk.c
 
$(LIB)(angle.o):	angle.c
 
$(LIB)(ascale.o):	ascale.c
 
$(LIB)(assign_name.o):	assign_name.c
 
$(LIB)(atomangle.o):	atomangle.c
 
$(LIB)(backbone_f.o):	backbone_f.c
 
$(LIB)(batch_job.o):	batch_job.c
 
$(LIB)(bin_search.o):	bin_search.c
 
$(LIB)(binschr.o):	binschr.f
 
$(LIB)(bondl.o):	bondl.c
 
$(LIB)(c_print_atom.o):	c_print_atom.c
 
$(LIB)(cartx2.o):	cartx2.f
 
$(LIB)(cclsngdst.o):	cclsngdst.c
 
$(LIB)(ccterxyz.o):	ccterxyz.c
 
$(LIB)(cgen.o):	cgen.c
 
$(LIB)(cgrdconsp.o):	cgrdconsp.c
 
$(LIB)(cgrdinit.o):	cgrdinit.c
 
$(LIB)(cgrdncons.o):	cgrdncons.c
 
$(LIB)(cgrdnxt.o):	cgrdnxt.c
 
$(LIB)(cgwrtinit.o):	cgwrtinit.c
 
$(LIB)(cgwrtnxt.o):	cgwrtnxt.c
 
$(LIB)(checknatc.o):	checknatc.c
 
$(LIB)(checkrtf.o):	checkrtf.c
 
$(LIB)(chkbbone.o):	chkbbone.c
 
$(LIB)(chkchnclsr.o):	chkchnclsr.c
 
$(LIB)(chkclsdst.o):	chkclsdst.c
 
$(LIB)(chkcntcts.o):	chkcntcts.c
 
$(LIB)(chkftres.o):	chkftres.c
 
$(LIB)(chmceil.o):	chmceil.c
 
$(LIB)(chnclsr_f.o):	chnclsr_f.c
 
$(LIB)(clear_atom.o):	clear_atom.c
 
$(LIB)(clear_atoms.o):	clear_atoms.c
 
$(LIB)(clear_side.o):	clear_side.c
 
$(LIB)(clschna.o):	clschna.f
 
$(LIB)(cnsbblist.o):	cnsbblist.c
 
$(LIB)(cnsfblist.o):	cnsfblist.c
 
$(LIB)(codes.o):	codes.c
 
$(LIB)(compspcgrd.o):	compspcgrd.f
 
$(LIB)(constatm.o):	constatm.c
 
$(LIB)(constclmp.o):	constclmp.c
 
$(LIB)(copyfarray.o):	copyfarray.c
 
$(LIB)(copyst.o):	copyst.c
 
$(LIB)(cpfbest.o):	cpfbest.c
 
$(LIB)(cpftbest.o):	cpftbest.c
 
$(LIB)(cpintobst.o):	cpintobst.c
 
$(LIB)(cprint.o):	cprint.f
 
$(LIB)(cprint1.o):	cprint1.c
 
$(LIB)(cpyprm.o):	cpyprm.c
 
$(LIB)(csearch.o):	csearch.c
 
$(LIB)(debug.o):	debug.c
 
$(LIB)(delatmfgrd.o):	delatmfgrd.c
 
$(LIB)(delatmfgrd1.o):	delatmfgrd1.f
 
$(LIB)(delatmsfgrd.o):	delatmsfgrd.c
 
$(LIB)(die.o):	die.c
 
$(LIB)(dispatch.o):	dispatch.c
 
$(LIB)(dot.o):	dot.c
 
$(LIB)(dststop.o):	dststop.c
 
$(LIB)(dump_grid.o):	dump_grid.c
 
$(LIB)(dwprangl.o):	dwprangl.f
 
$(LIB)(ecchkcont.o):	ecchkcont.f
 
$(LIB)(ecntrl.o):	ecntrl.c
 
$(LIB)(ecsetlims.o):	ecsetlims.f
 
$(LIB)(ecsetxyzsys.o):	ecsetxyzsys.f
 
$(LIB)(ephi.o):	ephi.c
 
$(LIB)(evaluate_f.o):	evaluate_f.c
 
$(LIB)(explcntcts.o):	explcntcts.c
 
$(LIB)(extnlims.o):	extnlims.c
 
$(LIB)(fill2.o):	fill2.c
 
$(LIB)(fill4.o):	fill4.c
 
$(LIB)(fill_grid.o):	fill_grid.c
 
$(LIB)(fill_grid1.o):	fill_grid1.f
 
$(LIB)(fillatm_c.o):	fillatm_c.c
 
$(LIB)(fillatm_ca.o):	fillatm_ca.c
 
$(LIB)(fillatm_cb.o):	fillatm_cb.c
 
$(LIB)(fillatm_h.o):	fillatm_h.c
 
$(LIB)(fillatm_ht1.o):	fillatm_ht1.c
 
$(LIB)(fillatm_n.o):	fillatm_n.c
 
$(LIB)(fillatm_o.o):	fillatm_o.c
 
$(LIB)(fillbbatms.o):	fillbbatms.c
 
$(LIB)(fillbndang.o):	fillbndang.c
 
$(LIB)(filllog.o):	filllog.f
 
$(LIB)(filspc.o):	filspc.c
 
$(LIB)(finbbone.o):	finbbone.c
 
$(LIB)(finchncls.o):	finchncls.c
 
$(LIB)(finclsatm.o):	finclsatm.c
 
$(LIB)(finres.o):	finres.c
 
$(LIB)(finsidchn.o):	finsidchn.c
 
$(LIB)(fixctero.o):	fixctero.c
 
$(LIB)(fixinitgrd.o):	fixinitgrd.c
 
$(LIB)(fixnterh.o):	fixnterh.c
 
$(LIB)(fixpdb.o):	fixpdb.c
 
$(LIB)(forprtatm.o):	forprtatm.f
 
$(LIB)(forscanf.o):	forscanf.c
 
$(LIB)(fortrace.o):	fortrace.f
 
$(LIB)(free_range.o):	free_range.c
 
$(LIB)(fresrk.o):	fresrk.c
 
$(LIB)(fwerf.o):	fwerf.c
 
$(LIB)(gammln.o):	gammln.f
 
$(LIB)(gammp.o):	gammp.f
 
$(LIB)(gammq.o):	gammq.f
 
$(LIB)(gcd.o):	gcd.c
 
$(LIB)(gcf.o):	gcf.f
 
$(LIB)(generate.o):	generate.c
 
$(LIB)(genh.o):	genh.c
 
$(LIB)(genic.o):	genic.c
 
$(LIB)(get_atnum.o):	get_atnum.c
 
$(LIB)(get_resnum.o):	get_resnum.c
 
$(LIB)(getclsangle.o):	getclsangle.c
 
$(LIB)(getclsbond.o):	getclsbond.c
 
$(LIB)(getdate.o):	getdate.c
 
$(LIB)(getgrdspc.o):	getgrdspc.c
 
$(LIB)(getparam.o):	getparam.c
 
$(LIB)(getparbond.o):	getparbond.f
 
$(LIB)(getparbond2.o):	getparbond2.c
 
$(LIB)(getres.o):	getres.c
 
$(LIB)(getseg.o):	getseg.c
 
$(LIB)(getstring.o):	getstring.c
 
$(LIB)(gettime.o):	gettime.c
 
$(LIB)(getuv.o):	getuv.f
 
$(LIB)(gser.o):	gser.f
 
$(LIB)(gtclschngom.o):	gtclschngom.c
 
$(LIB)(gtprangclsa.o):	gtprangclsa.c
 
$(LIB)(gtprangl.o):	gtprangl.f
 
$(LIB)(gtprangl2.o):	gtprangl2.c
 
$(LIB)(gtrbstatms.o):	gtrbstatms.c
 
$(LIB)(hadd.o):	hadd.c
 
$(LIB)(hbenergy.o):	hbenergy.f

FOBJS = \
abmerf.c abmerfc.c abmerfi.c binschr.c cartx2.c clschna.c compspcgrd.c  \
cprint.c delatmfgrd1.c dwprangl.c ecchkcont.c ecsetlims.c ecsetxyzsys.c  \
fill_grid1.c filllog.c forprtatm.c fortrace.c gammln.c gammp.c gammq.c  \
gcf.c getparbond.c getuv.c gser.c gtprangl.c hbenergy.c infrls.c  \
interpa.c maxr.c merfi.c mkprolrng.c nindx.c rnhbsrch.c rnnbsrch.c  \
schdls.c srchnat1.c srchnatm1.c stuphb1.c

csearch : $(FOBJS)
	


.f.o :
	f2c $(FOPT) $<
	mv $*.c $*.d
	sed 's/_(/(/' $*.d >$*.c
	rm $*.d
	cc -o $*.o -c $*.c
 


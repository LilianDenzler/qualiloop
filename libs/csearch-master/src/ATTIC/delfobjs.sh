FOBJS='abmerf abmerfc abmerfi binschr cartx2 clschna compspcgrd 
cprint delatmfgrd1 dwprangl ecchkcont ecsetlims ecsetxyzsys 
fill_grid1 filllog forprtatm fortrace gammln gammp gammq 
gcf getparbond getuv gser gtprangl hbenergy infrls 
interpa maxr merfi mkprolrng nindx rnhbsrch rnnbsrch 
schdls srchnat1 srchnatm1 stuphb1 main'

for fobj in $FOBJS
do
        rm $fobj.o
        rm $fobj.c
done





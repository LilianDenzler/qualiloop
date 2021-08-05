#ifndef __FZ_PRO__
#define __FZ_PRO__
 

/*
    Prototypes file for the Fz*.f modules
    Created by librar on 03-Nov-1992
*/

extern void ( * Fzadatm ) (
char   * AtomName   ,
FINT4  * NameSize1  ,
char   * Type       ,
FINT4  * TypeSize   ,
char   * Name       ,
FINT4  * NameSize2  ,
char   * Flags      ,
FINT4  * FlagSize   ,
FINT4  * AtomGroup  ,
FREAL8 * AtomCharge ,
FINT4  * error ) ;
 
extern void ( * Fzadbnd ) (
FINT4  * atom1 ,
FINT4  * atom2 ,
FINT4  * error ) ;
 
extern void ( * Fzadmol ) (
char   * Moltitle ,
FINT4  * MaxSize  ,
FINT4  * error ) ;
 
extern void ( * Fzadres ) (
char   * Name  ,
FINT4  * NameSize ,
char   * Number ,
FINT4  * NumberSize ,
FINT4  * error ) ;
 
extern void ( * Fzcrdst ) (
FINT4  * CoordSets ,
FINT4  * error ) ;
 
extern void ( * Fzdebug ) (
FINT4 * Flags ) ;
 
extern void ( * Fzet1rs ) (
FREAL8 * vscale ,
FINT4  * Error ) ;
 
extern void ( * Fzetcel ) (
FREAL8 * celldata ,
FINT4  * error ) ;
 
extern void ( * Fzetcnj ) (
FINT4  * formula , 
FREAL8 * tol , 
FREAL8 * sprec ,
FREAL8 * extbnd ,
FREAL8 * rfn , 
FREAL8 * der , 
FREAL8 * rms , 
FINT4  * Error ) ;
 
extern void ( * Fzetdyn ) (
FINT4  * isave   ,
FINT4  * restart ,
FINT4  * istep   ,
FINT4  * alg     ,
FINT4  * nstpav  ,
FREAL8 * timstp  ,
FREAL8 * demax   ,
FINT4  * error ) ;
 
extern void ( * Fzetene ) (
FINT4 * error ) ;
 
extern void ( * Fzeters ) (
FREAL8 * toten ,
FREAL8 * dtoten ,
FINT4  * error ) ;
 
extern void ( * Fzetffd ) (
FINT4  * morse  ,
FINT4  * cross  ,
FINT4  * oop    ,
FINT4  * inc1_4 ,
FINT4  * Ecut   ,
FREAL8 * cutdis ,
FREAL8 * swtdis ,
FINT4  * error ) ;
 
extern void ( * Fzetgea ) (
FINT4 * norder ,
FINT4 * ncorct ,
FINT4 * kclrh ,
FINT4 * error ) ;
 
extern void ( * Fzetidn ) (
FREAL8 * temp ,
FINT4  * step ,
FINT4  * error ) ;
 
extern void ( * Fzetkpr ) (
FREAL8 * press  ,
FREAL8 * time   ,
FREAL8 * compr  ,
FINT4  * kshape ,
FINT4  * kcnstr ,
FINT4  * katprs ,
FINT4  * katpfc ,
FINT4  * katscl ,
FINT4  * Error ) ;

 
extern void ( * Fzetktp ) (
FREAL8 * temp ,
FREAL8 * time ,
FINT4  * Error ) ;
 
extern void ( * Fzetmin ) (
FINT4 * save ,
FINT4 * steep ,
FINT4 * conj ,
FINT4 * qmin ,
FINT4 * newt ,
FINT4 * error ) ;
 
extern void ( * Fzetnbr ) (
FINT4  * nrneib ,
FREAL8 * cutoff ,
FREAL8 * xofadj ,
FINT4  * ReturnCode ) ;
 
extern void ( * Fzetnwt ) (
FREAL8 * accnr  ,
FREAL8 * exeig  ,
FREAL8 * stpnew ,
FINT4  * Error ) ;
 
extern void ( * Fzetqmn ) (
FINT4  * mode ,
FREAL8 * dfn ,
FREAL8 * hscale ,
FREAL8 * stepmin ,
FREAL8 * derstop ,
FREAL8 * rmsstop ,
FINT4  * error ) ;

 
extern void ( * Fzetrnd ) (
FREAL8 * temp ,
FREAL8 * seed ,
FINT4  * error ) ;
 
extern void ( * Fzetrsi ) (
char   * Restart ,
FINT4  * error ) ;
 
extern void ( * Fzetrso ) (
char   * Restart ,
FINT4  * error ) ;
 
extern void ( * Fzetstp ) (
FREAL8 * step    ,
FREAL8 * factinc ,
FREAL8 * factdec ,
FREAL8 * derstop ,
FREAL8 * rmsstop ,
FINT4  * error ) ;
 
extern void ( * Fzetsym ) (
FREAL8 * symmetry ,
FREAL8 * trans    ,
FINT4  * error ) ;
 
extern void ( * Fzettrs ) (
FREAL8 * temp   ,
FREAL8 * dtemp  ,
FINT4  * krsclv ,
FINT4  * kexact ,
FREAL8 * vscale ,
FINT4  * kindv  ,
FREAL8 * attemp ,
FREAL8 * atscal ,
FINT4  * error ) ;
 
extern void ( * Fzetxyz ) (
FINT4  * AtomUse ,
FREAL8 * X1      ,
FREAL8 * X2      ,
FREAL8 * X3      ,
FINT4  * error ) ;
 
extern void ( * Fzflcel ) (
FINT4  * nrbox  ,
FREAL8 * boxlim ,
FINT4  * ncellx ,
FINT4  * ncelly ,
FINT4  * ncellz ,
FINT4  * error ) ;
 
extern void ( * Fzfncon ) (
FINT4  * error ) ;
 
extern void ( * Fzfncrd ) (
FINT4  * error ) ;
 
extern void ( * Fzfndat ) (
FINT4  * error ) ;
 
extern void ( * Fzfnpar ) (
FINT4  * error ) ;
 
extern void ( * Fzfxmol ) (
FINT4  * molecule ,
FINT4  * error ) ;
 
extern void ( * Fzfxres ) (
FINT4  * residue ,
FINT4  * error ) ;
 
extern void ( * Fzfxtor ) (
FINT4  * Atom  ,
FREAL8 * angle ,
FREAL8 * force ) ;
 
extern void ( * Fzno14 ) (
FREAL8 * sf14lj ,
FREAL8 * sf14qq ,
FINT4  * error ) ;
 
extern void ( * Fznocrs ) (
FINT4  * error ) ;
 
extern void ( * Fznopln ) (
FINT4  * error ) ;
 
extern void ( * Fzpracc ) (
char   * Name     ,
FINT4  * NameSize ,
FINT4  * error ) ;
 
extern void ( * Fzprang ) (
char   * Name1    ,
char   * Name2    ,
char   * Name3    ,
FINT4  * NameSize ,
FREAL8 * Data     ,
FINT4  * error ) ;
 
extern void ( * Fzprbnd ) (
char   * Name1    ,
char   * Name2    ,
FINT4  * NameSize ,
FREAL8 * Data1    ,
FREAL8 * Data2    ,
FREAL8 * Data3    ,
FINT4  * error ) ;
 
extern void ( * Fzprchg ) (
char   * Name     ,
FINT4  * NameSize ,
FREAL8 * Charge   ,
FINT4  * error ) ;
 
extern void ( * Fzprhbd ) (
FREAL8 * Dist  ,
FINT4  * error ) ;
 
extern void ( * Fzprhyd ) (
char   * Name     ,
FINT4  * NameSize ,
FINT4  * error ) ;
 
extern void ( * Fzprmas ) (
char   * Name     ,
FINT4  * NameSize ,
FREAL8 * Mass     ,
FINT4  * error ) ;
 
extern void ( * Fzprnb2 ) (
char   * Name1    ,
char   * Name2    ,
FINT4  * NameSize ,
FREAL8 * Data     ,
FINT4  * error ) ;
 
extern void ( * Fzprnbd ) (
char   * Name     ,
FINT4  * NameSize ,
FREAL8 * Data     ,
FINT4  * error ) ;
 
extern void ( * Fzprpln ) (
char   * Name1    ,
char   * Name2    ,
char   * Name3    ,
char   * Name4    ,
FINT4  * NameSize ,
FREAL8 * Data     ,
FINT4  * error ) ;
 
extern void ( * Fzprtor ) (
char   * Name1    ,
char   * Name2    ,
char   * Name3    ,
char   * Name4    ,
FINT4  * NameSize ,
FREAL8 * Data     ,
FINT4  * error ) ;
 
extern void ( * Fzprtrf ) (
char   * Name1    ,
char   * Name2    ,
char   * Name3    ,
char   * Name4    ,
FINT4  * NameSize ,
FREAL8 * Data     ,
FINT4  * error ) ;
 
extern void ( * Fzprtrw ) (
char   * Name1    ,
char   * Name2    ,
FINT4  * NameSize ,
FREAL8 * Data     ,
FINT4  * error ) ;
 
extern void ( * Fzprtyp ) (
char   * Name ,
FINT4  * Size ,
FINT4  * icg  ,
FINT4  * iwg  ,
FINT4  * ibe  ,
FINT4  * ite  ,
FINT4  * ipe  ,
FINT4  * ioe  ,
FINT4  * inb  ,
FINT4  * iop  ,
FINT4  * ims  ,
FINT4  * error ) ;
 
extern void ( * Fzprxan ) (
char   * Name1    ,
char   * Name2    ,
char   * Name3    ,
char   * Name4    ,
FINT4  * NameSize ,
FREAL8 * Data     ,
FINT4  * error ) ;
 
extern void ( * Fzptcel ) (
FINT4  * nrbox  ,
FREAL8 * boxlim ,
FINT4  * ncellx ,
FINT4  * ncelly ,
FINT4  * ncellz ,
FINT4  * error ) ;
 
extern void ( * Fzstatn ) (
FINT4  * error ) ;


#endif

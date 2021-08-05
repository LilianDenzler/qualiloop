#ifndef __MISCFORT_PRO__
#define __MISCFORT_PRO__
 

#if defined (__FORT_UNDERSCORE__)
#define yadatm yadatm_
#endif
void     yadatm (
char   * AtomName ,
FINT4  * ATomNamelen ,
char   * Type ,
FINT4  * Typelen ,
char   * Name ,
FINT4  * Namelen ,
char   * Flags ,
FINT4  * Flagslen ,
FINT4  * Group ,
FREAL8 * Charge ,
FINT4  * ErrorCode ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yadbnd yadbnd_
#endif
void     yadbnd  (
FINT4  * Atom1 ,
FINT4  * Atom2 ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yadmol yadmol_
#endif
void     yadmol (
char   * Line ,
FINT4  * Size ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yadres yadres_
#endif
void     yadres (
char   * Name ,
FINT4  * SizeName ,
char   * Number ,
FINT4  * SizeNumb ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define ycrdst ycrdst_
#endif
void     ycrdst (
FINT4  * CoordSets ,
FINT4  * error
) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yetcel yetcel_
#endif
void     yetcel (
FREAL8 * cellvar ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yetsym yetsym_
#endif
void     yetsym (
FREAL8 * symmetry ,
FREAL8 * translation ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yetxyz yetxyz_
#endif
void     yetxyz  (
FINT4  * Atom ,
FREAL8 * x ,
FREAL8 * y ,
FREAL8 * z ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yflcel yflcel_
#endif
void     yflcel (
FINT4  * nrbox ,
FREAL8 * boxlim ,
FINT4  * ncellx ,
FINT4  * ncelly ,
FINT4  * ncellz ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yfndat yfndat_
#endif
void     yfndat (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yfxmol yfxmol_
#endif
void    yfxmol (
FINT4 * Molecule ,
FINT4 * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yfxres yfxres_
#endif
void    yfxres (
FINT4 * Residue ,
FINT4 * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yfxtor yfxtor_
#endif
void     yfxtor (
FINT4  * atoms ,
FREAL8 * angle ,
FREAL8 * force ,
FINT4  * error ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yno14 yno14_
#endif
void     yno14  (
FREAL8 * sf14lj ,
FREAL8 * sf14qq ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define ynocrs ynocrs_
#endif
void     ynocrs (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define ynopln ynopln_
#endif
void     ynopln (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define yptcel yptcel_
#endif
void     yptcel (
FINT4  * nrbox ,
FREAL8 * boxlim ,
FINT4  * ncellx ,
FINT4  * ncelly ,
FINT4  * ncellz ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define ystatn ystatn_
#endif
void     ystatn (
FINT4  * Code ) ;

#if defined (__FORT_UNDERSCORE__)
#define zadatm zadatm_
#endif
void     zadatm (
char   * AtomName ,
FINT4  * ATomNamelen ,
char   * Type ,
FINT4  * Typelen ,
char   * Name ,
FINT4  * Namelen ,
char   * Flags ,
FINT4  * Flagslen ,
FINT4  * Group ,
FREAL8 * Charge ,
FINT4  * ErrorCode ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zadbnd zadbnd_
#endif
void     zadbnd  (
FINT4  * Atom1 ,
FINT4  * Atom2 ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zadmol zadmol_
#endif
void     zadmol (
char   * Line ,
FINT4  * Size ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zadres zadres_
#endif
void     zadres (
char   * Name ,
FINT4  * SizeName ,
char   * Number ,
FINT4  * SizeNumb ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zdebug zdebug_
#endif
void     zdebug (
FINT4  * kdebug ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zet1rs zet1rs_
#endif
void     zet1rs (
FREAL8 * vscale ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetcel zetcel_
#endif
void     zetcel (
FREAL8 * cellvar ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetcnj zetcnj_
#endif
void     zetcnj (
FINT4  * formula ,
FREAL8 * tol ,
FREAL8 * sprec ,
FREAL8 * extbnd ,
FREAL8 * rfn  ,
FREAL8 * derstp ,
FREAL8 * rmsstp ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetdyn zetdyn_
#endif
void     zetdyn (
FINT4  * nsave ,
FINT4  * steps ,
FINT4  * constPT ,
FINT4  * algorithm ,
FINT4  * Vrescale ,
FINT4  * nstpav ,
FREAL8 * timstep ,
FREAL8 * demax ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetene zetene_
#endif
void     zetene (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zeters zeters_
#endif
void     zeters (
FREAL8 * toten ,
FREAL8 * dtoten ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetffd zetffd_
#endif
void     zetffd (
FINT4  * JMORSE ,
FINT4  * JCROSS ,
FINT4  * JOPLNG ,
FINT4  * JEXC   ,
FINT4  * JCUT   ,
FREAL8 * CUTDSI ,
FREAL8 * SWTDSI ,
FINT4  * IERROR ) ;

 
#if defined (__FORT_UNDERSCORE__)
#define zetgea zetgea_
#endif
void     zetgea (
FINT4  * order ,
FINT4  * correct ,
FINT4  * clear_h ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetidn zetidn_
#endif
void zetidn (
FINT4  * TEMPIN ,
FREAL8 * JSTEP  ,
FINT4  * IERROR ) ;

 
#if defined (__FORT_UNDERSCORE__)
#define zetkpr zetkpr_
#endif
void     zetkpr (
FREAL8 * press ,
FREAL8 * time ,
FREAL8 * compress ,
FINT4  * kshape ,
FINT4  * kcnstr ,
FINT4  * katprs ,
FINT4  * katpfc ,
FINT4  * katscl ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetktp zetktp_
#endif
void     zetktp (
FREAL8 * temp ,
FREAL8 * time ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetmin zetmin_
#endif
void     zetmin (
FINT4  * isave ,
FINT4  * isteep ,
FINT4  * iconj ,
FINT4  * iqmin ,
FINT4  * inewt ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetnbr zetnbr_
#endif
void     zetnbr (
FINT4  * nrneib ,
FREAL8 * cutoff ,
FREAL8 * xofadj ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetnwt zetnwt_
#endif
void     zetnwt (
FREAL8 * accnr ,
FREAL8 * exeig ,
FREAL8 * stpnew ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetqmn zetqmn_
#endif
void     zetqmn (
FINT4  * mode ,
FREAL8 * dfn ,
FREAL8 * hscale ,
FREAL8 * exstm ,
FREAL8 * derstp ,
FREAL8 * rmsstp ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetrnd zetrnd_
#endif
void     zetrnd (
FREAL8 * temp ,
FREAL8 * seed ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetrsi zetrsi_
#endif
void     zetrsi (
FINT4  * JRSTRT ,
FINT4  * IERROR ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetrso zetrso_
#endif
void     zetrso (
FINT4  * JRSTRT ,
FINT4  * IERROR ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetstp zetstp_
#endif
void     zetstp (
FREAL8 * step0 ,
FREAL8 * fctinc ,
FREAL8 * fctdec ,
FREAL8 * derstp ,
FREAL8 * rmsstp ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetsym zetsym_
#endif
void     zetsym (
FREAL8 * symmetry ,
FREAL8 * translation ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zettrs zettrs_
#endif
void     zettrs (
FREAL8 * temp ,
FREAL8 * dtemp ,
FINT4  * kisclv ,
FINT4  * kexact ,
FREAL8 * vscale ,
FINT4  * kindv ,
FREAL8 * attemp ,
FREAL8 * atscal ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zetxyz zetxyz_
#endif
void     zetxyz  (
FINT4  * Atom ,
FREAL8 * x ,
FREAL8 * y ,
FREAL8 * z ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zflcel zflcel_
#endif
void     zflcel (
FINT4  * nrbox ,
FREAL8 * boxlim ,
FINT4  * ncellx ,
FINT4  * ncelly ,
FINT4  * ncellz ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zfncon zfncon_
#endif
void     zfncon (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zfncrd zfncrd_
#endif
void     zfncrd (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zfndat zfndat_
#endif
void     zfndat (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zfnpar zfnpar_
#endif
void     zfnpar (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zfxmol zfxmol_
#endif
void    zfxmol (
FINT4 * Molecule ,
FINT4 * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zfxres zfxres_
#endif
void    zfxres (
FINT4 * Residue ,
FINT4 * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zfxtor zfxtor_
#endif
void     zfxtor (
FINT4  * atoms ,
FREAL8 * angle ,
FREAL8 * force ,
FINT4  * error ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zno14 zno14_
#endif
void     zno14  (
FREAL8 * sf14lj ,
FREAL8 * sf14qq ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define znocrs znocrs_
#endif
void     znocrs (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define znopln znopln_
#endif
void     znopln (
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zpracc zpracc_
#endif
void     zpracc (
char   * Name ,
FINT4  * NameSize ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprang zprang_
#endif
void     zprang (
char   * Name1 ,
char   * Name2 ,
char   * Name3 ,
FINT4  * NameSize ,
FREAL8 * Pot ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprbnd zprbnd_
#endif
void     zprbnd (
char   * Name1 ,
char   * Name2 ,
FINT4  * NameSize ,
FREAL8 * Pot1 ,
FREAL8 * Pot2 ,
FREAL8 * Pot3 ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprchg zprchg_
#endif
void     zprchg (
char   * Name ,
FINT4  * NameSize ,
FREAL8 * Charge ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprhbd zprhbd_
#endif
void     zprhbd (
FREAL8 * Dist ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprhyd zprhyd_
#endif
void     zprhyd (
char   * Name ,
FINT4  * NameSize ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprmas zprmas_
#endif
void     zprmas (
char   * Name ,
FINT4  * NameSize ,
FREAL8 * Mass ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprnb2 zprnb2_
#endif
void     zprnb2 (
char   * Name1 ,
char   * Name2 ,
FINT4  * NameSize ,
FREAL8 * Pot ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprnbd zprnbd_
#endif
void     zprnbd (
char   * Name ,
FINT4  * NameSize ,
FREAL8 * Pot ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprpln zprpln_
#endif
void     zprpln (
char   * Name1 ,
char   * Name2 ,
char   * Name3 ,
char   * Name4 ,
FINT4  * NameSize ,
FREAL8 * Pot ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprtor zprtor_
#endif
void     zprtor (
char   * Name1 ,
char   * Name2 ,
char   * Name3 ,
char   * Name4 ,
FINT4  * NameSize ,
FREAL8 * Pot ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprtrf zprtrf_
#endif
void     zprtrf (
char   * INAME1 ,
char   * INAME2 ,
char   * INAME3 ,
char   * INAME4 ,
FINT4  * ILEN   ,
FREAL8 * POTS   ,
FINT4  * IERROR ) ;

 
#if defined (__FORT_UNDERSCORE__)
#define zprtrw zprtrw_
#endif
void     zprtrw (
char   * Name1 ,
char   * Name2 ,
FINT4  * NameSize ,
FREAL8 * Pot ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprtyp zprtyp_
#endif
void     zprtyp (
char   * Name ,
FINT4  * NameSize ,
FINT4  * icg ,
FINT4  * iwg ,
FINT4  * ibe ,
FINT4  * ite ,
FINT4  * ipe ,
FINT4  * ioe ,
FINT4  * inb ,
FINT4  * iop ,
FINT4  * ims ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zprxan zprxan_
#endif
void     zprxan (
char   * Name1 ,
char   * Name2 ,
char   * Name3 ,
char   * Name4 ,
FINT4  * NameSize ,
FREAL8 * Pot ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zptcel zptcel_
#endif
void     zptcel (
FINT4  * nrbox ,
FREAL8 * boxlim ,
FINT4  * ncellx ,
FINT4  * ncelly ,
FINT4  * ncellz ,
FINT4  * Code ) ;
 
#if defined (__FORT_UNDERSCORE__)
#define zstatn zstatn_
#endif
void     zstatn (
FINT4  * Code ) ;

#endif

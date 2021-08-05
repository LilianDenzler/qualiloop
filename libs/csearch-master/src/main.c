#include "ProtoTypes.h"
#include "CongenProto.h"

#ifdef DEFINE_HERE
/* Common Block Declarations */
struct {
    integer natoms, nres, nsegs, nbonds, nangs, nptors, nitors, nhbs, nnbs, 
	    ndonat, naccat, nbpar, napar, nptpar, nitpar, nhbpar, natyps, 
	    nbauto;
} values_;

#define values_1 values_

struct {
    real grid2, cutnb2, cuthb2, epsilon;
    logical cons_die__;
    real cos_cuthba__;
    integer ioff[100], parm_nop__, aamapp, glymapp, promapp, naamap, nglymap, 
	    npromap, proconsp, proconsphip, eproconsp, nprocons;
    char ctitle[800]	/* was [80][10] */;
    integer nctitl, glymapu, alamapu, promapu, proconsu, stunit;
    real glyemax, alaemax, proemax, eringpro;
    logical ignoreevdw;
    integer maxleaf, restart_stlen__, restart_st__;
    logical save_coor__;
    integer nconsp, consp, savex, savey, savez;
    logical sidehits_opt__;
    real maxdt_def__;
} cg_;

#define cg_1 cg_

struct {
    integer dof_headp__;
} dof_head_;

#define dof_head__1 dof_head__

struct {
    integer dof_tailp__;
} dof_tail_;

#define dof_tail__1 dof_tail__

struct {
    integer confnum;
} confnum_;

#define confnum_1 confnum_

struct {
    integer leafnum;
} leafnum_;

#define leafnum_1 leafnum_

struct {
    integer top_level_envp__;
} top_level_env_;

#define top_level_env__1 top_level_env__

struct {
    integer dummy_sidehitsp__;
} dummy_sidehits_;

#define dummy_sidehits__1 dummy_sidehits__

struct {
    real xcart[6150], ycart[6150], zcart[6150], xwork[6150], ywork[6150], 
	    zwork[6150];
} coords_;

#define coords_1 coords_

struct {
    integer dbg_cgen__, dbg_clschn__, dbg_alloc__, dbg_allstk__, dbg_allhp__;
} dbg_;

#define dbg_1 dbg_

struct {
    real eqbdis[150], bndcon[150], eqang[350], angcon[350], torphs[75], 
	    tormlt[75], torcon[75], eqitan[55], impcon[55], vdwr12[1640], 
	    vdwr6[1640], hbr12[250], hbr10[250], atmpol[100], atneff[100], 
	    vdwrad[100];
    shortint atflag[100];
    integer bndkey[150], angkey[350], impkey[55], torkey[75], nbkey[1640], 
	    hbkey[250], nbcut, dielec, nbflag;
} engpar_;

#define engpar_1 engpar_

struct {
    integer grid_space_grid__, ngridx, ngridy, ngridz;
    real xmn, ymn, zmn, xmx, ymx, zmx, spgridsz, recipgrid;
    integer grid_excluded__, grid_cntnbx__, grid_ingrid__, maxnbx, 
	    grid_nbxa__, freehd, freecls, grid_nexthd__, grid_clshd__, 
	    grid_clstl__, grid_nextcls__, grid_clsatm__, grid_donp__, 
	    grid_accp__, grid_resbya__, grid_qside__, grid_radius__;
    real maxradius;
} grid_;

#define grid_1 grid_

struct {
    shortint ibndp[6250], iangp[9150], itorp[3600], iimpp[3250], ihbp[4100];
} getpar_;

#define getpar_1 getpar_

struct {
    shortint hbdon[4100], hbacc[4100];
    real hbcut, hbacut;
} hbonds_;

#define hbonds_1 hbonds_

struct {
    shortint atbnd1[6250], atbnd2[6250], atang1[9150], atang2[9150], atang3[
	    9150], attor1[3600], attor2[3600], attor3[3600], attor4[3600], 
	    atimp1[3250], atimp2[3250], atimp3[3250], atimp4[3250], hbacpt[
	    1200], hbaan1[1200], hbaan2[1200], hbdonr[1200], hbdhyd[1200], 
	    hbdan1[1200], hbdan2[1200], lstatm[1051], resndx[1050], nbexcl[
	    16150], qmove[6150], atcode[6150];
    integer lstexc[6150], segndx[210]	/* was [10][21] */;
    real segid[20], resid[1050], resnme[1050], atmnme[6150], atchrg[6150], 
	    atmass[6150];
} pstruct_;

#define pstruct_1 pstruct_

struct {
    shortint bndat1[4000]	/* was [100][40] */, bndat2[4000]	/* 
	    was [100][40] */, angat1[6000]	/* was [150][40] */, angat2[
	    6000]	/* was [150][40] */, angat3[6000]	/* was [150][
	    40] */, torat1[1400]	/* was [35][40] */, torat2[1400]	
	    /* was [35][40] */, torat3[1400]	/* was [35][40] */, torat4[
	    1400]	/* was [35][40] */, impat1[880]	/* was [22][40] */, 
	    impat2[880]	/* was [22][40] */, impat3[880]	/* was [22][40] */, 
	    impat4[880]	/* was [22][40] */, exclnb[5600]	/* was [140][
	    40] */, nexcld[2800]	/* was [70][40] */, acpthb[320]	/* 
	    was [8][40] */, aan1hb[320]	/* was [8][40] */, aan2hb[320]	/* 
	    was [8][40] */, donrhb[320]	/* was [8][40] */, dan1hb[320]	/* 
	    was [8][40] */, dan2hb[320]	/* was [8][40] */, dhydhb[320]	/* 
	    was [8][40] */, fstgrp[400]	/* was [10][40] */, lstgrp[400]	/* 
	    was [10][40] */, acindx[2800]	/* was [70][40] */, nparam[
	    400]	/* was [10][40] */, nparus[400]	/* was [10][40] */, 
	    nparmx[10];
    integer nreses;
    real atmmas[100], atmchg[70], b1bld[2800]	/* was [70][40] */, b2bld[
	    2800]	/* was [70][40] */, a1bld[2800]	/* was [70][40] */, 
	    a2bld[2800]	/* was [70][40] */, torbld[2800]	/* was [70][
	    40] */, bldat1[2800]	/* was [70][40] */, bldat2[2800]	
	    /* was [70][40] */, bldat3[2800]	/* was [70][40] */, bldat4[
	    2800]	/* was [70][40] */, namres[40], acodes[100], namatm[
	    2800]	/* was [70][40] */, grpnam[400]	/* was [10][40] */;
} restop_;

#define restop_1 restop_

struct {
    real sc_bond_bld__[200], sc_angle_bld__[200], sc_tors_bld__[200], 
	    sc_offset__[200];
    integer nscatm, nsclmp, nscres, sc_code_bld__[200], sc_atom_part__[76], 
	    sc_clump_part__[31], sc_special__[30], sc_symmetry__[75], 
	    sc_ante1_bld__[200], sc_ante2_bld__[200], sc_ante3_bld__[200], 
	    sc_resname__[30], sc_free_atom__[75], sc_atom_bld__[200];
} sidetop_;

#define sidetop_1 sidetop_
#endif 

int main(int argc, char **argv)
{
    csearch(argc, argv);
    return(0);
}


#ifndef __CG_H__
#define __CG_H__

#define first 1
#define independent 2
#define all 3
#define combination 4
#define iterative 5
 
#define eo_energy 1
#define eo_rms 2
 
#define eval_code_rms 1
#define eval_code_energy 2
#define eval_code_user 3

#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif
 
EXTERN struct dof *dof_head,*dof_tail;
 
EXTERN int confnum,leafnum;
 
EXTERN jmp_buf *top_level_env;
 
EXTERN int *dummy_sidehits;
 
EXTERN struct {
    float grid2,cutnb2,cuthb2,epsilon;
    logical cons_die;
    float cos_cuthba;
    int ioff[100];
/*    short int *parm_no;  ACRM 20.02.06 */
    short int parm_no[100][100];
    struct pepmap *aamap,*glymap,*promap;
    int naamap,nglymap,npromap;
    float (*procons)[9],
	  *proconsphi,
	  (*eprocons)[3];
    int nprocons;
    char ctitle[10][80];
    int nctitl;
    int glymapu,alamapu,promapu,proconsu,stunit;
    float glyemax,alaemax,proemax,eringpro;
    logical ignore_evdw;
    int maxleaf;
    int restart_stlen;
    char *restart_st;
    logical save_coor;
    int nconsp;
    int *consp;
    float *savex,*savey,*savez;
    logical sidehits_opt;
    float maxdt_def;
    } cg;

#endif

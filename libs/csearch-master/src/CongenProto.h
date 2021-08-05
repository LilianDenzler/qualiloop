/* 08.04.94  Added prototype for clear_sc()     By: ACRM
*/

/* This is the new CONGEN prototype/headers file for CR 318 */
/* Created by: John Woods, Oxford Molecular Ltd, 11/08/92 */
/* It contains prototypes for all CONGEN routines, plus some
other headers previously defined in the file standards.h */
/* NB: All changes to Congen routine prototypes and prototypes
for new Congen routines *must* be appended to this file */

/* OMUPD DW 19/08/92 Add underscores to certain structures etc. if required */

#ifdef __FORT_UNDERSCORE__
#define cg cg_
#define hbonds hbonds_
#define nbnd nbnd_
#define values values_
#define dbg dbg_
#define pstruct pstruct_
#define restop restop_
#define dof_head dof_head_
#define top_level_env top_level_env_
#define engpar engpar_
#define getpar getpar_
#define coords coords_
#define dof_tail dof_tail_
#define confnum confnum_
#define leafnum leafnum_
#define dummy_sidehits dummy_sidehits_
#define sidetop sidetop_
#define nindx nindx_
#define grid grid_
#define stuphb1 stuphb1_
#define mkprolrng mkprolrng_
#define srchnatm1 srchnatm1_
#define infrls infrls_
#endif

#include <setjmp.h>
#include "params.h"
#include "cgentypes.h"
#include "cg.h"
#include "cgen.h"
#include "cgmacros.h"
#include "getpar.h"
#include "conglob.h"
#include "coords.h"
#include "ctitla.h"
#include "dbg.h"
#include "forscanf.h"
#include "grid.h"
#include "hbonds.h"
#include "macros.h"
#include "numpointers.h"
#include "engpar.h"
#include "parse.h"
#include "parse_file.h"
#include "pdb.h"
#include "pstruct.h"
#include "restop.h"
#include "restopio.h"
#include "sdef.h"
#include "sidetop.h"

/* Defines for C routines called from Fortran follow */

#ifdef __FORT_UNDERSCORE__
#define ascale ascale_
#define cprint1 cprint1_
#define dot dot_
#define die die_
#define fwerf fwerf_
#define generate generate_
#define initfgener initfgener_
#define init_debug init_debug_
#define initfcgen initfcgen_
#define pscgncmnd pscgncmnd_
#define prtdbgvars prtdbgvars_
#define padspace padspace_
#define print_atom1 print_atom1_
#define phia phia_
#define parse parse_
#define regrid regrid_
#define ReadRTF ReadRTF_
#define ReadParams ReadParams_
#define ReadCoords ReadCoords_
#define rbest_f rbest_f_
#define setup_parser setup_parser_
#define SetupSeq SetupSeq_
#define setdefs setdefs_
#define schnatm schnatm_
#define sidechain_f sidechain_f_
#define status_f status_f_
#define wrtcoordsf wrtcoordsf_
#define vlen vlen_
#define termcgen termcgen_

/* OMUPD DW 18/08/92 Added missing functions which need _ appending */

#define fill2 fill2_
#define fill4 fill4_
#define adatmtgrd adatmtgrd_
#define fixinitgrd fixinitgrd_
#define cpyprm cpyprm_
#define trace trace_
#define chmceil chmceil_
#define c_print_atom c_print_atom_
#define getseg getseg_

#endif


/* Now follow C prototypes for all C and Fortran routines */


#if defined (__FORT_UNDERSCORE__)
#define abmerf abmerf_
#endif
/* No prototype for abmerf - please supply one */
#if defined (__FORT_UNDERSCORE__)
#define abmerfc abmerfc_
#endif
/* No prototype for abmerfc - please supply one */
#if defined (__FORT_UNDERSCORE__)
#define abmerfi abmerfi_
#endif
/* No prototype for abmerfi - please supply one */
int adafls ( 
int *head , 
int *tail , 
int *p , 
int *new , 
int *next
 ) 
;
void adatmstgrd ( 
struct atom **atoms
 ) 
;
void adatmtgrd ( 
int ind
 ) 
;
void adccangs ( 
int **cangle, 
float cmaxdt[] , 
int *ianglep , 
struct chain_closure_d *desc
 ) 
;
int AddH ( 
PDB *hlist , 
PDB **position , 
int igtype_m ) 
;
char *alloc ( 
int n
 ) 
;
void allstk  ( 
int numwrd
 ) 
;
void angle  ( 
short *iarray , 
short *jarray , 
short *karray , 
logical *tarray , 
logical *itar , 
int *nangle , 
single *aarray , 
single *x , 
single *y , 
single *z
 ) 
;
float angleposneg ( float opp , 
                  float adj ) 
;
float atomangle ( float xi , 
                float yi , 
                float zi , 
                float xj , 
                float yj , 
                float zj , 
                float xk , 
                float yk , 
                float zk ) 
;
void RotateVec ( float *x ,  
               float *y ,  
               float *z ,  
               float angle ,  
               int   axis ) 
;
float phi ( float x1 , 
          float y1 , 
          float z1 , 
          float x2 , 
          float y2 , 
          float z2 , 
          float x3 , 
          float y3 , 
          float z3 , 
          float x4 , 
          float y4 , 
          float z4 ) 
;
float trueangle ( float opp , 
                float adj ) 
;

void ascale ( 
float *u , 
float *scl , 
float *v , 
int *dim
 ) 
;
void assign_name  ( 
char  type[maxat][4] , 
short *ibase , 
int   *ires , 
int   *index_var , 
char  alpha_var[]
 ) 
;
void backbone_f  ( 
struct dof *dofp , 
struct backbone_d *desc
 ) 
;
int batch_job ( 
void
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define binschr binschr_
#endif
/* No prototype for binschr - please supply one */
void bondl  ( 
short *iarray , 
short *jarray , 
short *karray , 
logical *tarray , 
short *nbonds , 
single *x , 
single *y , 
single *z , 
single *barray
 ) 
;
int c_nindx ( 
int numb , 
int *narray , 
int nlen
 ) 
;
void c_print_atom ( 
char string[] , 
int  atnum
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define cartx2 cartx2_
#endif
/* No prototype for cartx2 - please supply one */
void cclsngdst  ( 
struct backbone_d *desc , 
int *catom , 
int natom , 
logical stretch_all , 
int **cangle, 
float cmaxdt[] , 
int nangle
 ) 
;
int ccterxyz  ( 
PDB *p ,       
PDB *ca_p ,    
PDB *c_p ,     
PDB *o_p     
 ) 
;
void cgen  ( 
void
 ) 
;
void cgrdconsp  ( 
FILE *fp , 
int  *consp , 
int  *ncons
 ) 
;
void cgrdinit  ( 
FILE *fp , 
int  *ok
 ) 
;
void cgrdncons  ( 
FILE *fp , 
int  *ncons
 ) 
;
void cgrdnxt  ( 
FILE *fp , 
int *consp , 
int ncons , 
float *tote , 
int *eof
 ) 
;
void cgwrtinit  ( 
int *consp , 
int ncons
 ) 
;
void cgwrtnxt  ( 
int *consp , 
int ncons , 
float tote
 ) 
;
void CheckNATC ( 
int *natc
 ) 
;
void CheckRTF ( 
void
 ) 
;
void chkbbone ( 
int startres , 
int lastres , 
logical *okp
 ) 
;
void chkchnclsr ( 
int startres , 
int lastres , 
logical *okp
 ) 
;
logical chkclsdst ( 
struct backbone_d *desc
 ) 
;
logical chkcntcts ( 
struct atom **atomsp , 
float *maxevdw , 
int *sidehits
 ) 
;
void chkftres ( 
struct backbone_d *desc
 ) 
;
int chmceil ( 
float *r
 ) 
;
void chnclsr_f ( 
struct dof *dofp , 
struct chain_closure_d *desc
 ) 
;
void clear_atom ( 
struct atom *p
 ) 
;
void clear_atoms ( 
char chain[] , 
char res1[] , 
char res2[]
 ) 
;

void clear_side(char *res1, char *res2);
void clear_sc(char *chain, char *res1, char *res2);

#if defined (__FORT_UNDERSCORE__)
#define clschna clschna_
#endif
/* No prototype for clschna - please supply one */
void cnsbblist ( 
struct backbone_d *desc , 
int *catom , 
int *natomp
 ) 
;
void cnsfblist ( 
struct backbone_d *desc , 
int *catom , 
int *natomp
 ) 
;
void codes ( 
int lpsf , 
int lhb
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define compspcgrd compspcgrd_
#endif
/* No prototype for compspcgrd - please supply one */
void constatm ( 
struct atom *atp , 
float torsion
 ) 
;
void constclmp ( 
struct clump *clp , 
float phi
 ) 
;
void CopyFArray ( 
float *a , 
float *cop , 
int   n
 ) 
;
void copyst ( 
char st1[] , 
int *st1siz , 
int *st1len , 
char st2[] , 
int *st2len
 ) 
;
void cpfbest ( 
struct sideres *srp , 
int ind
 ) 
;
void cpftbest ( 
struct sideres *srp , 
int fromind , 
int toind
 ) 
;
void cpintobst ( 
struct sideres *srp , 
int ind
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define cprint cprint_
#endif
/* No prototype for cprint - please supply one */
void cprint1 ( 
char string[] , 
int *len
 ) 
;
void cpyprm ( 
float *parm , 
short *itc , 
short *iac , 
float *atoma , 
int *natom
 ) 
;
void csearch ( 
int argc , 
char **argv
 ) 
;
void debug ( 
void
 ) 
;
void delatmfgrd ( 
int ind
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define delatmfgrd1 delatmfgrd1_
#endif
/* No prototype for delatmfgrd1 - please supply one */
void delatmsfgrd ( 
struct atom **atoms
 ) 
;
void die ( 
void
 ) 
;
void dispatch ( 
struct dof *dofp
 ) 
;
float dot ( 
float *u , 
float *v , 
int *dim
 ) 
;
void dststop ( 
void
 ) 
;
void dump_grid ( 
void
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define dwprangl dwprangl_
#endif
/* No prototype for dwprangl - please supply one */
#if defined (__FORT_UNDERSCORE__)
#define ecchkcont ecchkcont_
#endif
/* No prototype for ecchkcont - please supply one */
int ecntrl ( 
void
 ) ;
#if defined (__FORT_UNDERSCORE__)
#define ecsetlims ecsetlims_
#endif
/* No prototype for ecsetlims - please supply one */
#if defined (__FORT_UNDERSCORE__)
#define ecsetxyzsys ecsetxyzsys_
#endif
/* No prototype for ecsetxyzsys - please supply one */
double ephi ( 
int itype , 
short *ip , 
short *jp , 
short *kp , 
short *lp , 
short *icp , 
int nphi , 
float *cpc , 
float *cpd , 
float *cpb , 
float *x , 
float *y , 
float *z , 
short *imove , 
int iprint , 
int imaxp
 ) ;

void evaluate_f ( 
struct dof *dofp , 
struct evaluate_d *desc) ;

void explcntcts ( 
struct sideres *srp , 
struct clump *clp , 
struct range_list **rangepp , 
int *sidehits
 ) 
;
void extnlims ( 
void
 ) 
;
void fill2 ( 
short *a , 
int *n , 
short *value
 ) 
;
void fill4 ( 
long *a , 
int *n , 
long *value
 ) 
;
void fill_grid ( 
void
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define fill_grid1 fill_grid1_
#endif

void fillatm_c ( 
struct atom **atp , 
int resno , 
logical forward
 ) 
;
void fillatm_ca ( 
struct atom **atp , 
int resno , 
logical forward
 ) 
;
void fillatm_cb ( 
struct atom **atp , 
int resno , 
logical forward
 ) 
;
void fillatm_h ( 
struct atom **atp , 
int resno , 
logical forward
 ) 
;
void fillatm_ht1 ( 
struct atom **atp , 
int resno , 
logical forward
 ) 
;
void fillatm_n ( 
struct atom **atp , 
int resno , 
logical forward
 ) 
;
void fillatm_o ( 
struct atom **atp , 
int resno , 
logical forward
 ) 
;
void fillbbatms ( 
struct atom **atoms , 
int resno , 
logical forward , 
logical nter , 
struct atom **at_n , 
struct atom **at_h , 
struct atom **at_ca , 
struct atom **at_cb , 
struct atom **at_cg , 
struct atom **at_cd , 
struct atom **at_c , 
struct atom **at_o , 
struct atom **at_ht1
 ) 
;
void fillbndang ( 
struct atom *atp
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define filllog filllog_
#endif
/* No prototype for filllog - please supply one */
void filspc ( 
char st[] , 
int  *stmax , 
int *stlen
 ) 
;
void finbbone ( 
struct dof *dofp , 
struct backbone_d *desc
 ) 
;
void finchncls ( 
struct dof *dofp , 
struct chain_closure_d *dp
 ) 
;
void finclsatm ( 
struct backbone_d *desc
 ) 
;
void finres ( 
float totchg , 
int nbxn[] , 
char nbxnam[mxnber][4] , 
int ntabl , 
char tblkey[MXTABL][4] , 
int *p_sawtbl , 
int tblrep[MXTABL]
 ) 
;
void finsidchn ( 
struct dof *dofp , 
struct sidechain_d *desc , 
int *res , 
float *sgrid , 
float *maxevdw , 
logical *vavoid , 
logical *clump_symmetry
 ) 
;
void FixCterO ( 
PDB *pdb
 ) 
;
void fixinitgrd ( 
void
 ) 
;
void FixNterH ( 
PDB *pdb
 ) 
;
void fixpdb ( 
FILE *hadd_fp , 
PDB *pdb
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define forprtatm forprtatm_
#endif
/* No prototype for forprtatm - please supply one */
int forscanf ( 
FILE *fp , 
char string[] , 
FLIST *list
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define fortrace fortrace_
#endif
/* No prototype for fortrace - please supply one */
void free_range ( 
struct range_list *rp
 ) 
;
int frestk ( 
int numwrd
 ) 
;
void fwerf ( float *xx ,  float *xtemp ) 
;
#if defined (__FORT_UNDERSCORE__)
#define gammln gammln_
#endif
/* No prototype for gammln - please supply one */
#if defined (__FORT_UNDERSCORE__)
#define gammp gammp_
#endif
/* No prototype for gammp - please supply one */
#if defined (__FORT_UNDERSCORE__)
#define gammq gammq_
#endif
/* No prototype for gammq - please supply one */
int gcd ( 
int in , 
int jn
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define gcf gcf_
#endif
/* No prototype for gcf - please supply one */
void generate ( 
char segname[] , 
int nbxmod
 ) 
;
int GenH ( 
PDB *pdb , 
unsigned short *err_flag , 
char grname[8]
 ) 
;
void genic ( 
int istart , 
int istop
 ) 
;
int get_atnum ( 
char chain[] , 
char res[] , 
char atm[]
 ) 
;
int get_resnum ( 
char chain[] , 
char res[]
 ) 
;
float getclsangle ( 
int i , 
int j , 
int k
 ) 
;
float getclsbond ( 
int i , 
int j
 ) 
;
void GetDate ( 
char day[]
 ) 
;
void getgrdspc  ( 
void
 ) 
;
int GetParam ( 
char  command[] , 
float *value , 
int   *nletters
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define getparbond getparbond_
#endif
/* No prototype for getparbond - please supply one */
float getparbond2 ( 
int i , 
int j
 ) 
;
int getres ( 
short atom , 
short *ibase , 
int nres
 ) 
;
int getseg ( 
int *ires , 
int nictot[maxseg+1][10] , 
int *nseg
 ) 
;
int GetString ( 
char command[] , 
char strparam[]
 ) 
;
void GetTime ( 
char time[]
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define getuv getuv_
#endif
/* No prototype for getuv - please supply one */
#if defined (__FORT_UNDERSCORE__)
#define gser gser_
#endif
/* No prototype for gser - please supply one */
void gtclschngom ( 
struct chain_closure_d *dp
 ) 
;
float gtprangclsa ( 
int it , 
int jt , 
int kt , 
logical stretch_all , 
int **cangle, 
float cmaxdt[] , 
int nangle
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define gtprangl gtprangl_
#endif
/* No prototype for gtprangl - please supply one */
float gtprangl2 ( 
int i , 
int j , 
int k
 ) 
;
void gtrbstatms ( 
struct dof *dofp , 
struct rbest_d *desc , 
logical *okp
 ) 
;
int hadd ( 
FILE *fp , 
PDB  *pdb
 ) 
;
#if defined (__FORT_UNDERSCORE__)
#define hbenergy hbenergy_
#endif
/* No prototype for hbenergy - please supply one */
void indexf ( 
int n , 
float *arrin , 
int   *indx
 ) 
;
void indexi ( 
int n , 
int *arrin , 
int *indx
 ) 
;
void init_debug ( 
void
 ) 
;
void initfcgen  ( 
void
 ) 
;
void initfgener ( 
void
 ) 
;
void isctrnglst ( 
struct range_list *in1 , 
struct range_list *in2 , 
struct range_list *sout
 ) 
;
char *killspcs ( 
char string[]
 ) 
;
void ljust ( 
char string[]
 ) 
;
void ljustpad ( 
char string[]
 ) 
;
int LookupName ( 
char name[] , 
int  *ok , 
int  ntabl , 
char tblkey[MXTABL][4] , 
int  tblrep[MXTABL] , 
int  *sawtbl
 ) 
;
PDB *makeh ( 
int igtype_m , 
float gr_m , 
float alpha_m , 
float beta_m , 
unsigned short firstres
 ) 
;
void makinb ( 
int mode
 ) 
;
void mapgrdtorng ( 
float lbound , 
float hbound , 
float grid , 
struct range_list *rangep , 
float **philistp , 
int *nphip
 ) 
;
int match ( 
char comstring[] , 
char string2[] , 
int  *nletters
 ) 
;
int matom ( 
int nres , 
char atom[]
 ) 
;
float maxcoor ( 
float *a , 
int   n
 ) 
;
float maxf2 ( 
float x1 , 
float x2
 ) 
;
int maxi2 ( 
int x1 , 
int x2
 ) 
;
float mincoor ( 
float *a , 
int   n
 ) 
;
float minf2 ( 
float x1 , 
float x2
 ) 
;
float minf3 ( 
float x1 , 
float x2 , 
float x3
 ) 
;
int mini2 ( 
int x1 , 
int x2
 ) 
;
float mygetparbond ( 
int *ib , 
int *jb , 
short int parm_no[100][100] , 
short int *iac , 
int   *kcb , 
int   *ncb , 
char  type[maxat][4] , 
float *cbb
 ) 
;
void mynextwd ( 
char *instring[] , 
char word[]
 ) 
;
void onethr ( 
char one , 
char three[]
 ) 
;
void pack2 ( 
unsigned long int *i4 , 
unsigned short int i2 , 
unsigned short int j2
 ) 
;
void pad ( 
char string[] , 
int num
 ) 
;
void padspace ( 
char string[] , 
int  length
 ) 
;
void padterm ( 
char string[] , 
int size
 ) 
;
int parse ( 
char *comline , 
int nkeys , 
KeyWd *keywords , 
float *floatparam , 
char *strparam[]
 ) 
;
void patexh ( 
int ires , 
int iatom , 
int iphi , 
int imphi , 
int ndons
 ) 
;
void phia ( 
short *iarray , 
short *jarray , 
short *karray , 
short *larray , 
single *parray , 
int *ndih , 
single *x , 
single *y , 
single *z
 ) 
;
void prdie ( 
char string[]
 ) 
;
void print_atom ( 
struct atom *atp
 ) 
;
void print_atom1 ( 
char st[] , 
int iat
 ) 
;
void print_atoms ( 
struct atom **atpp
 ) 
;
void procang ( 
FILE *fp , 
short IndxTab[100][100]
 ) 
;
void ProcessBonds ( 
FILE *fp , 
short IndxTab[100][100]
 ) 
;
void prochbnd ( 
FILE *fp , 
short IndxTab[100][100]
 ) 
;
void procimps ( 
FILE *fp , 
short IndxTab[100][100]
 ) 
;
void procnbnds ( 
FILE *fp , 
short IndxTab[100][100] , 
int *ncn
 ) 
;
void proctors ( 
FILE *fp , 
short IndxTab[100][100]
 ) 
;
void propbbones ( 
struct dof **p , 
struct backbone_d *dp , 
int startres , 
int lastres
 ) 
;
void prtcgncmnd ( 
void
 ) 
;
void prtdbgvars ( 
void
 ) 
;
void prtglobopt ( 
void
 ) 
;
void prtgrdsz ( 
void
 ) 
;
void prtrnglst ( 
struct range_list *rp
 ) 
;
void prtsidchn ( 
struct sidechain_d *dp
 ) 
;
int pscgncmnd ( 
void
 ) 
;
void rbest_f ( 
struct dof *dofp , 
struct rbest_d *desc
 ) 
;
void rdemap ( 
int unit , 
int *nmap , 
struct pepmap **map , 
float emax , 
char name[]
 ) 
;
void rdpepmaps ( 
void
 ) 
;
void rdprocons ( 
int unit , 
float  ( **procons ) [9] , 
float **proconsphi , 
float  ( **eprocons ) [3] , 
int *nprocons , 
float eringpro
 ) 
;
void rdprtitle ( 
FILE *fp
 ) 
;
void ReadCoords ( 
FILE *fp , 
FILE *pgp_fp , 
int showcoords
 ) 
;
int ReadParams ( 
FILE *fp
 ) 
;
void ReadPDB ( 
FILE *fp , 
PDB *pdb , 
int *natom
 ) 
;
int ReadPGP ( 
FILE *fp
 ) 
;
int ReadPIR ( 
FILE *fp , 
int  maxaa , 
char *seqs[]
 ) 
;
void ReadRTF ( 
FILE *unit , 
int qprint , 
int *natc
 ) 
;
void regrid ( 
int ind
 ) 
;
void RTFAcceptor ( 
int *p_sawtbl , 
int ntabl , 
char tblkey[MXTABL][4] , 
int *tblrep
 ) 
;
void RTFAngle ( 
int *p_sawtbl , 
int ntabl , 
char tblkey[MXTABL][4] , 
int *tblrep
 ) 
;
void RTFAtom ( 
int *nbxn , 
int sawtbl , 
int *nrx , 
char nbxnam[mxnber][4] , 
int *natc
 ) 
;
void RTFBond ( 
int *p_sawtbl , 
int ntabl , 
char tblkey[MXTABL][4] , 
int tblrep[MXTABL]
 ) 
;
void RTFBuild ( 
void
 ) 
;
void RTFDeclare ( 
char tblkey[MXTABL][4] , 
int *tblrep , 
int *ntabl
 ) 
;
void RTFDonor ( 
int *p_sawtbl , 
int ntabl , 
char tblkey[MXTABL][4] , 
int *tblrep
 ) 
;
void RTFGroup ( 
int *p_sawtbl , 
int ntabl , 
char tblkey[MXTABL][4] , 
int *tblrep
 ) 
;
void RTFImproper ( 
int *p_sawtbl , 
int ntabl , 
char tblkey[MXTABL][4] , 
int *tblrep
 ) 
;
void RTFPrint ( 
int *prnflg
 ) 
;
void RTFResidue ( 
float *p_totchg , 
int *nrx , 
int *p_sawtbl , 
int *nbxn , 
char nbxnam[mxnber][4] , 
int ntabl , 
char tblkey[MXTABL][4] , 
int tblrep[MXTABL]
 ) 
;
int RTFTopology ( 
void
 ) 
;
void RTFTorsion ( 
int *p_sawtbl , 
int ntabl , 
char tblkey[MXTABL][4] , 
int *tblrep
 ) 
;
void RTFType ( 
int *natc
 ) 
;
void schnatm ( 
int ind , 
int mode , 
float *maxevdw , 
logical ignore_evdw , 
int *sidehits , 
logical *impact , 
float *eel , 
float *evdw , 
float *ehb
 ) 
;
void scnemapln ( 
char *line , 
float *omega , 
float *phi , 
float *psi , 
float *e1 , 
float *e2
 ) 
;
void selsidetor ( 
struct sideres *srp , 
struct resphi **philist , 
short ip[] , 
short jp[] , 
short kp[] , 
short lp[] , 
short icp[] , 
int nphi
 ) 
;
void setavdclp ( 
struct sideres *srp , 
struct clump *clp
 ) 
;
void setchncls ( 
float *maxdt , 
float *maxg , 
float *maxevdw , 
logical *cistrans , 
int *startres , 
logical *ok
 ) 
;
void setdefs ( 
void
 ) 
;
void setrcntcts ( 
struct atom *atp , 
float maxevdw , 
short parm_no[100][100]
 ) 
;
void  setrestitr ( 
void
 ) 
;
void setsdbound ( 
struct sideres *srp , 
struct clump *clp , 
float *lbound , 
float *hbound , 
float *grid
 ) 
;
void setsddofrs ( 
int startres , 
int lastres , 
int *res , 
float def_sgrid , 
float def_maxevdw , 
logical def_vavoid , 
logical def_clump_symmetry , 
int **presp , 
float **psgridp , 
float **pmaxevdwp , 
logical **pvavoidp , 
logical **pclump_symmetryp
 ) 
;
void setsidchn ( 
int prev_start , 
int prev_last , 
struct sidechain_d *dp , 
int *res , 
float *sgrid , 
float *maxevdw , 
logical *vavoid , 
logical *clump_symmetry , 
logical *okp
 ) 
;
void setsidphis ( 
struct sideres *srp , 
struct clump *clp , 
float **philistp , 
int *nphip , 
int *sidehits
 ) 
;
void setup_consp ( 
int **conspp , 
int *nconspp , 
struct dof *last_dof
 ) 
;
void setup_hbond ( 
void
 ) 
;
void setup_imove ( 
struct evaluate_d *desc
 ) 
;
void setup_nbond ( 
void
 ) 
;
void setup_parmno ( 
short parm_no[100][100]
 ) 
;
void setup_parser ( 
void
 ) 
;
void setup_status ( 
struct status_d *desc , 
logical *okp
 ) 
;
void SetupSeq ( 
int nchain , 
char *sequence
 ) 
;
void shufcombs ( 
struct dof *dofp , 
struct sidechain_d *desc , 
struct sideres **srpp
 ) 
;
int  sidchnall ( 
struct dof *dofp , 
struct sidechain_d *desc , 
struct sideres **srpp , 
struct clump **clpp
 ) 
;
void sidchnbst ( 
struct sidechain_d *desc , 
int maxconf , 
struct sideres *srp , 
struct clump **clpp , 
logical never_fail
 ) 
;
logical sidchncomb ( 
struct dof *dofp , 
struct sidechain_d *desc
 ) 
;
double sidchnengy ( 
struct sideres *srp
 ) 
;
int  sidchnfst ( 
struct dof *dofp , 
struct sidechain_d *desc , 
logical dispatch_flag , 
struct sideres **srpp , 
struct clump **clpp
 ) 
;
logical sidchnindi ( 
struct dof *dofp , 
struct sidechain_d *desc
 ) 
;
void sidchnitr ( 
struct dof *dofp , 
struct sidechain_d *desc
 ) 
;
double sidchnrms ( 
struct sideres *srp
 ) 
;
void sidechain_f ( 
struct dof *dofp , 
struct sidechain_d *desc
 ) 
;
int sign ( 
int i
 ) 
;
void skiptitle ( 
FILE *fp
 ) 
;
void sort_i ( 
int *arr , 
int n
 ) 
;
void sort_ui ( 
unsigned long int *arr , 
int n
 ) 
;
void sorti_perm ( 
int *array , 
int numb , 
int *idx
 ) 
;
int srchint ( 
int *array , 
int numval , 
int val
 ) 
;
int srchwd ( 
char words[][4] , 
int  numwrd , 
char word[]
 ) 
;
int srwdbd ( 
char a[][4] , 
int *low , 
int *high , 
char wd[]
 ) 
;
void statprn ( 
int count , 
int *a , 
int n
 ) 
;
void status_f ( 
struct dof *dofp , 
struct status_d *desc
 ) 
;
void stread ( 
FILE *stunit
 ) 
;
void struppr ( 
char string1[] , 
char string2[]
 ) 
;
void termcgen ( 
void
 ) 
;

void terminate(char *string);

void trace ( 
void
 ) 
;
void trimic ( 
void
 ) 
;
int typeinres ( 
int ires , 
char iupac[]
 ) 
;
void unpack2 ( 
unsigned long int i4 , 
unsigned short int *i2 , 
unsigned short int *j2
 ) 
;
void updtnictot ( 
int mode , 
int prflag
 ) 
;
int user_cgeval ( 
void
 ) 
;
float vlen ( 
float *v , 
int   *dim
 ) 
;
void write_status ( 
int count , 
int *iters , 
int niter , 
int unit , 
int *write_count , 
char fst[] , 
int fstmax , 
int fstlen
 ) 
;
void WritePDB ( 
FILE *fp , 
PDB *pdb
 ) 
;
void wrtcoordsf ( 
struct dof *dofp , 
struct write_coordinates_d *desc
 ) 
;
void wrtpdbrec ( 
FILE *fp , 
PDB *pdb
 ) 
;
void xbyind ( 
int asize , 
char a[] , 
char b[] , 
int aind[] , 
int nind
 ) 
;

/* ACRM 17.02.06 added prototypes
 */
void compspcgrd(int *ind, int *ix, int *iy, int *iz, logical *oob);
void mkprolrng(float *x, float *y, float *z,
               int *ind_c, int *ind_n, int *ind_ca, int *ind_cb, 
               int *ind_cg, int *ind_gd, float *phi, float (*procons)[],
               float *proconsphi, int *nprocons);
void cartx2(float *x, float *y, float *z, int *i1, int *i2, int *i3, int *i4,
            float *bond, float *theta, float *phi);
void clschna(float *x, float *y, float *z,
             int (*atmind)[], float (*bond)[], float (*angle)[],
             float *newx, float *newy, float *newz, float *omega,
             int *iter, logical *entzus, float *pomega, float *maxdt,
             float (*nangle)[], float *maxg);
void delatmfgrd1(int *ind, short int *space_grid,
                 logical *ingrid, logical *nexthd,
                 int *clshd, int *clstl, int *nextcls, int *clsatm);
void ecsetxyzsys(int *ante, real *ori, real *ex, real *ey, real *ez);
void ecsetlims(float *oc, float *ecx, float *avoidxcenter, float *avoidextent,
               float *avoidrmax, float *avoidvdwrmax, int *startx, int *starty,
               int *startz, int *lastx, int *lasty, int *lastz, float *rmax);

void ecchkcont(int *iatm, int *atomno, float *ori, 
               float *ecx, float *ecy, float *ecz, 
               float *avoidx, float *avoidr, float *avoidphi,
               float *rcontact, float *startphi,
               float *lastphi, logical *impact, int *sidehits,
               int *resbya, int *cntnbx, short int *nbxa, logical *qside);
/*
void fill_grid1(short int *space_grid, logical *excluded,
                int *cntnbx, logical *ingrid, short int *nbxa,
                int *nexthd, int *clshd, int *clstl, int *nextcls,
                int *clsatm, int *resbya, float *radius);
*/
int nindx(int *number, int *narray, int *nlen);
void abmpad(char string[], int num);
void  dwprangl(int *it, int *jt, int *kt, short int parm_no[100][100],
               short int *iac, int *kct, int *nct, int *natc, char type[maxat][4],
               float *ctb, float *retval);
int bin_search(int numb, int *narray, int nlen);
void  infrls(int *head, int *tail, int *next, int *n);
void srchnatm1(int *ind, short int *space_grid,
               int *search_mode, float *maxevdw,
               logical *ignore_evdw, int *cntnbx, short int *nbxa,
               logical *excluded,
               short int parm_no[100][100],
               int *clshd, int *nextcls, int *clsatm, int *donp,
               int *accp, logical *qside, int *sidehits, int *resbya, logical *impact,
               double *eel, double *evdw, double *ehb, logical *outofbound);
void stuphb1(int *donp, int *accp);


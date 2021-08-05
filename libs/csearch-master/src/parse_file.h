#ifndef __PARSE_FILE_H__
#define __PARSE_FILE_H__

#define MAXCOMLEN 200

/************************************************************************
Set up global keywords for the command parser
*************************************************************************/

#define NCOMM 37                   /* Number of keywords                */
#define MAXNUMPARAM 10             /* Max number of numeric parameters  */
#define MAXSTRPARAM 10             /* Max number of string parameters   */
#define MAXSTRLEN   80             /* Max length of a string param      */
 
#ifdef MAIN
KeyWd   g_keys[NCOMM];             /* Array of keywords                 */
char   *g_strparam[MAXSTRPARAM];   /* Array of returned strings         */
float   g_numparam[MAXNUMPARAM];   /* Array of returned floats          */
#else
extern KeyWd   g_keys[NCOMM];             /* Array of keywords          */
extern char   *g_strparam[MAXSTRPARAM];   /* Array of returned strings  */
extern float   g_numparam[MAXNUMPARAM];   /* Array of returned floats   */
#endif

/************************************************************************/

/* Setup defines for parser */

#define do_chain    0
#define do_forward  1
#define do_reverse  2
#define do_side     3
#define do_write    4
#define do_cgen     5
#define do_loopnum  6
#define do_restop   7
#define do_params   8
#define do_topology 9
#define do_coords   10
#define do_sidetop  11
#define do_alamap   12
#define do_glymap   13
#define do_promap   14
#define do_procons  15
#define do_status   16
#define do_restart  17
#define do_loops    18
#define do_nbond    19
#define do_glyemax  20
#define do_alaemax  21
#define do_proemax  22
#define do_eringpro 23
#define do_debug    24
#define do_eps      25
#define do_cutnb    26
#define do_evaluate 27
#define do_cuthb    28
#define do_cutha    29
#define do_sequence 30
#define do_pgp      31
#define do_echo     32
#define do_end      33
#define do_clear    34
#define do_nofix    35
#define do_sclear   36
 
/************************************************************************/
/* Global file pointers */
#ifdef MAIN
FILE *in         = NULL,   /* The input command file     */
     *out        = NULL,   /* Output listing file        */
     *glymap_fp  = NULL,   /* Glycine map                */
     *alamap_fp  = NULL,   /* Alanine map                */
     *promap_fp  = NULL,   /* Proline map                */
     *procons_fp = NULL,   /* Proline constructor        */
     *status_fp  = NULL,   /* Status file                */
     *sidetop_fp = NULL,   /* Sidechain topology         */
     *loops_fp   = NULL,   /* Database loops for restart */
     *write_fp   = NULL,   /* Output conformations file  */
     *restop_fp  = NULL,   /* Residue topology file      */
     *parm_fp    = NULL,   /* Parameter file             */
     *pgp_fp     = NULL,   /* Proton generation params   */
     *pdb_fp     = NULL,   /* Coordinate PDB             */
     *seq_fp     = NULL;   /* Sequence file              */
#else
extern FILE *in,           /* The input command file     */
            *out,          /* Output listing file        */
            *glymap_fp,    /* Glycine map                */
            *alamap_fp,    /* Alanine map                */
            *promap_fp,    /* Proline map                */
            *procons_fp,   /* Proline constructor        */
            *status_fp,    /* Status file                */
            *sidetop_fp,   /* Sidechain topology         */
            *loops_fp,     /* Database loops for restart */
            *write_fp,     /* Output conformations file  */
            *restop_fp,    /* Residue topology file      */
            *parm_fp,      /* Parameter file             */
            *pgp_fp,       /* Proton generation params   */
            *pdb_fp,       /* Coordinate PDB             */
            *seq_fp;       /* Sequence file              */
#endif

#endif

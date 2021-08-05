#include "ProtoTypes.h"
#include "CongenProto.h"
 
void setup_parser(void)
{
   int i;
 
   /* Initialise string array */
   for(i=0; i<MAXSTRPARAM; i++)
   {
      g_strparam[i] = (char *)alloc(MAXSTRLEN * sizeof(char));
   }
 
   /* Construct keywords */
   /* For CGEN options */
   MAKEKEY(g_keys[do_cgen],    "CGEN",        NUMBER,0);
   MAKEKEY(g_keys[do_chain],   "CHAIN",       STRING,3);
   MAKEKEY(g_keys[do_forward], "FORWARD",     STRING,6);
   MAKEKEY(g_keys[do_reverse], "REVERSE",     STRING,6);
   MAKEKEY(g_keys[do_side],    "SIDECHAIN",   STRING,6);
   MAKEKEY(g_keys[do_write],   "WRITE",       STRING,1);
   MAKEKEY(g_keys[do_loops],   "LOOPS",       STRING,1);
   MAKEKEY(g_keys[do_restart], "RESTART",     STRING,1);
   MAKEKEY(g_keys[do_end],     "END",         NUMBER,0);
   /* For global options */
   MAKEKEY(g_keys[do_loopnum], "LOOPNUM",     NUMBER,1);
   MAKEKEY(g_keys[do_restop],  "RESTOP",      STRING,1);
   MAKEKEY(g_keys[do_params],  "PARAMS",      STRING,1);
   MAKEKEY(g_keys[do_topology],"TOPOLOGY",    STRING,1);
   MAKEKEY(g_keys[do_coords],  "COORDS",      STRING,1);
   MAKEKEY(g_keys[do_sidetop], "SIDETOP",     STRING,1);
   MAKEKEY(g_keys[do_alamap],  "ALAMAP",      STRING,1);
   MAKEKEY(g_keys[do_glymap],  "GLYMAP",      STRING,1);
   MAKEKEY(g_keys[do_promap],  "PROMAP",      STRING,1);
   MAKEKEY(g_keys[do_procons], "PROCONS",     STRING,1);
   MAKEKEY(g_keys[do_status],  "STATUS",      STRING,1);
   MAKEKEY(g_keys[do_nbond],   "NBOND",       NUMBER,10);
   MAKEKEY(g_keys[do_glyemax], "GLYEMAX",     NUMBER,1);
   MAKEKEY(g_keys[do_alaemax], "ALAEMAX",     NUMBER,1);
   MAKEKEY(g_keys[do_proemax], "PROEMAX",     NUMBER,1);
   MAKEKEY(g_keys[do_eringpro],"ERINGPRO",    NUMBER,1);
   MAKEKEY(g_keys[do_debug],   "DEBUG",       STRING,2);
   MAKEKEY(g_keys[do_eps],     "EPS",         NUMBER,1);
   MAKEKEY(g_keys[do_cutnb],   "CUTNB",       NUMBER,1);
   MAKEKEY(g_keys[do_cuthb],   "CUTHB",       NUMBER,1);
   MAKEKEY(g_keys[do_cutha],   "CUTHA",       NUMBER,1);
   MAKEKEY(g_keys[do_sequence],"SEQUENCE",    STRING,1);
   MAKEKEY(g_keys[do_pgp],     "PGP",         STRING,1);
   MAKEKEY(g_keys[do_echo],    "ECHO",        NUMBER,0);
   MAKEKEY(g_keys[do_clear],   "CLEAR",       STRING,3);
   MAKEKEY(g_keys[do_nofix],   "NOHADD",      NUMBER,0);
   /* OMUPD rkw 23/11/93 Added SCLEAR to clear sidechains */
   MAKEKEY(g_keys[do_sclear],  "SCLEAR",      STRING,3);
#ifdef EVALUATE
   MAKEKEY(g_keys[do_evaluate],"EVALUATE",    NUMBER,0);
#endif
}
 

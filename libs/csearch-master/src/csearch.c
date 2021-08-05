#define MAIN 1
#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
#include "cgparse.h"

/* Main program */
 
int wrnlev,info_level;
int fixpdbfile = TRUE;
FILE *fp_suc;

void csearch(
int argc,
char **argv
)
{
   UBYTE eof = FALSE;
   char  comline[MAXCOMLEN+1];
   int   key, i,
         nseq = 0,
         showcoords = 0;
   char  *seqs[8],
         chain[8][4];
#ifdef F2C
   char  arg[80];
   int   one = 1,
         two = 2;
#endif

/* OMUPD rkw 03/11/92 Added a header */
printf("\n+++ MODULE: CONGEN +++\n");
printf(  "======================\n");


/* ACRM 07.01.94 Modified to read and write from stdin/stdout */
/*
//   if(argc < 3)
//   {
//      printf("Usage: csearch <file.inp> <file.out>\n");
//      exit(0);
//   }
*/

/* ACRM 25.03.94 Now reads files from command line if specified
   stdin/stdout otherwise. Handles f2c as special case.
*/
 
   /* Set the default input and output files */
   in = stdin;
   out = stdout;

#ifdef F2C
   getarg_(&one, arg, 80);
   if(isalpha(arg[0]))
   {
      int i;
      for(i=0;i<80 && arg[i] != ' ';i++) ;
      if(i<80) arg[i] = '\0';

      if((in=fopen(arg,"r"))==NULL)
      {
         printf("unable to open input file: `%s'\n",arg);
         exit(1);
      }
   }

   getarg_(&two, arg, 80);
   if(isalpha(arg[0]))
   {
      int i;
      for(i=0;i<80 && arg[i] != ' ';i++) ;
      if(i<80) arg[i] = '\0';

      if((out=fopen(arg,"w"))==NULL)
      {
         printf("Unable to open output file: `%s'\n",arg);
         exit(1);
      }
   }
#else
   if(argc > 1)
   {
      if((in=fopen(argv[1],"r"))==NULL)
      {
         printf("Unable to open input file: %s\n",argv[1]);
         exit(1);
      }
   }

   if(argc > 2)
   {
      if((out=fopen(argv[2],"w"))==NULL)
      {
         printf("Unable to open output file: %s\n",argv[2]);
         exit(1);
      }
   }
#endif
 
   /* Setup the command parser */
   setup_parser();

   setdefs();

   /* ACRM 18.01.06 - Allocate memory for restart string! */
   cg.restart_st = alloc(80*sizeof(char));
   cg.restart_st[0] = '\0';

   /* Initialise the debugging variables */
   init_debug();

   /* Setup default chain labels */
   strcpy(chain[0],"L");
   strcpy(chain[1],"H");
   strcpy(chain[2],"A");
   strcpy(chain[3],"B");
   strcpy(chain[4],"C");
   strcpy(chain[5],"D");
   strcpy(chain[6],"E");
   strcpy(chain[7],"F");
 
   /* Sit in a loop which reads the input file and parses commands */
   while(!eof && !feof(in))
   {
      if(!fgets(comline,MAXCOMLEN,in)) break;
      terminate(comline);
 
      /* Echo to output */
      fprintf(out,"CSEARCH> %s\n",comline);
 
      key = parse(comline,NCOMM,g_keys,g_numparam,g_strparam);
      switch(key)
      {
      case PARSE_ERRC:
         fprintf(out,"Unknown command: %s\n",comline);
         break;
      case PARSE_ERRP:
         fprintf(out,"Incorrect parameters in command: %s\n",comline);
         break;
      case PARSE_COMMENT:
         /* Do nothing */
         break;
      case do_cgen:
         cgen();     /* The main conf. search routine */
         break;
      case do_glymap:
         if((glymap_fp = fopen(g_strparam[0],"r"))==NULL)
         {
            fprintf(out,"Warning: Unable to open glymap: %s\n",g_strparam[0]);
         }
         break;
      case do_alamap:
         if((alamap_fp = fopen(g_strparam[0],"r"))==NULL)
         {
            fprintf(out,"Warning: Unable to open alamap: %s\n",g_strparam[0]);
         }
         break;
      case do_promap:
         if((promap_fp = fopen(g_strparam[0],"r"))==NULL)
         {
            fprintf(out,"Warning: Unable to open promap: %s\n",g_strparam[0]);
         }
         break;
      case do_procons:
         if((procons_fp = fopen(g_strparam[0],"r"))==NULL)
         {
            fprintf(out,"Warning: Unable to open procons: %s\n",g_strparam[0]);
         }
         break;
      case do_sidetop:
         if((sidetop_fp = fopen(g_strparam[0],"r"))==NULL)
         {
            fprintf(out,"Warning: Unable to open sidetop: %s\n",g_strparam[0]);
         }
         break;
      case do_status:
         if((status_fp = fopen(g_strparam[0],"r"))==NULL)
         {
            fprintf(out,"Warning: Unable to open status file: %s\n",
                    g_strparam[0]);
         }
         break;
      case do_alaemax:
         cg.alaemax = g_numparam[0];
         break;
      case do_glyemax:
         cg.glyemax = g_numparam[0];
         break;
      case do_proemax:
         cg.proemax = g_numparam[0];
         break;
      case do_eringpro:
         cg.eringpro = g_numparam[0];
         break;
      case do_restart:
         strcpy(cg.restart_st,g_strparam[0]);
         cg.restart_stlen = strlen(cg.restart_st);
         break;
      case do_debug:
         debug();
         break;
      case do_eps:
         engpar.dielec = g_numparam[0];
         break;
      case do_cutnb:
         engpar.nbcut = g_numparam[0];
         break;
      case do_cuthb:
         hbonds.hbcut = g_numparam[0];
         break;
      case do_cutha:
         hbonds.hbacut = g_numparam[0];
         break;
      case do_restop:
         if((restop_fp=fopen(g_strparam[0],"r"))==NULL)
            fprintf(out,"Warning: Unable to open restop: %s\n",g_strparam[0]);
         else
            ReadRTF(restop_fp,1,&values.natyps);
         break;
      case do_params:
         if((parm_fp=fopen(g_strparam[0],"r"))==NULL)
            fprintf(out,"Warning: Unable to open parameters: %s\n",
                         g_strparam[0]);
         else
            ReadParams(parm_fp);
         break;
      case do_sequence:
         if((seq_fp=fopen(g_strparam[0],"r"))==NULL)
            fprintf(out,"Warning: Unable to open sequence: %s\n",g_strparam[0]);
         else if(!restop_fp)
            fprintf(out,"Warning: Residue topology must be specified first\n");
         else if(!parm_fp)
            fprintf(out,"Warning: Parameter file must be specified first\n");
         else
         {
            nseq = ReadPIR(seq_fp,maxres,seqs);
     /* For each chain, set up the sequence and generate PSTRUCT segment */
            for(i=0; i<nseq; i++)
            {
               SetupSeq(i+1,seqs[i]);
               generate(chain[i],2); /* 2 is the nbxmod update mode */
            }
         }
         /* Now set up the codes arrays */
         codes(1,1);
         break;
      case do_pgp:
         if((pgp_fp=fopen(g_strparam[0],"r"))==NULL)
            fprintf(out,"Warning: Unable to open PGP: %s\n",g_strparam[0]);
         break;
      case do_coords:
         if((pdb_fp=fopen(g_strparam[0],"r"))==NULL)
            fprintf(out,"Warning: Unable to open PDB: %s\n",g_strparam[0]);
         else if(!pgp_fp)
            fprintf(out,"Warning: PGP must be specified first\n");
         else if(!seq_fp)
            fprintf(out,"Warning: Sequence must be specified first\n");
         else if(!restop_fp)
            fprintf(out,"Warning: Residue topology must be specified first\n");
         else
            ReadCoords(pdb_fp,pgp_fp,showcoords);
         break;
      case do_nbond:
         printf("!!!   NBOND keyword not yet handled   !!!\n");
         break;
      case do_echo:
         showcoords = 1;
         break;
      case do_nofix:
         fixpdbfile = FALSE;
         break;
      case do_clear:
         if(!pdb_fp) fprintf(out,"Warning: PDB must be specified first\n");
         padspace(g_strparam[0],4);
         padspace(g_strparam[1],4);
         padspace(g_strparam[2],4);
         clear_atoms(g_strparam[0],g_strparam[1],g_strparam[2]);
         break;
      case do_sclear:    /* Added ACRM 28.03.94                          */
         if (!pdb_fp)
         {
            fprintf(out,"Warning: PDB must be specified first\n");
         }
         padspace(g_strparam[0],4);
         padspace(g_strparam[1],4);
         padspace(g_strparam[2],4);
         clear_sc(g_strparam[0],g_strparam[1],g_strparam[2]);
         break;
      default:
         fprintf(out,"main() parser: Case not handled\n");
         break;
      }
   }
/* OMUPD rkw 02/11/92 Output SUCCESS.AbM file for error checking */

if((fp_suc = fopen("SUCCESS.AbM","w")) == NULL)
   {
   printf("ERROR - could not open file SUCCESS.AbM\n");
   exit(1);
   }
else
   {
   if (confnum == 1)
      fprintf(fp_suc,"CONGEN finished normally, producing %d conformation\n",
           confnum);
   else
      fprintf(fp_suc,"CONGEN finished normally, producing %d conformations\n",
           confnum);
  
   fclose(fp_suc);
   }

}
 

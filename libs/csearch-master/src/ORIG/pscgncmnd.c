#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* The main command parser in main() has hit a CGEN command, so
   switches control to this routine */
 
int pscgncmnd(
void
)
{
   int stmax = mxcbuf;
   char st[mxcbuf+1]; /* Extra char for null */
   char buf[81];
   int ind,stlen,i;
   struct dof *p = NULL;
   logical parse_ok,ok;
   int startres,lastres;
   char comline[MAXCOMLEN+1];
   int   eof = FALSE;
   int   key;
 
   parse_ok = true;
 
   cg.restart_st = alloc(80*sizeof(char));
   cg.restart_st[0] = '\0';
 
   if(!sidetop_fp)
      prdie("You must specify the sidechain topology file before CGEN\n");
 
   stread(sidetop_fp);    /* Has to come before sidechain processing. */
 
   while(!eof && !feof(in))
   {
      if(!fgets(comline,MAXCOMLEN,in)) break;
      terminate(comline);
 
      /* Echo to output */
      fprintf(out,"CGEN>    %s\n",comline);
 
      key = parse(comline,NCOMM,g_keys,g_numparam,g_strparam);
 
      switch(key)
      {
      case PARSE_ERRC:
         fprintf(out,"Unknown command: %s\n",comline);
         continue;
         break;
      case PARSE_ERRP:
         fprintf(out,"Incorrect parameters in command: %s\n",comline);
         continue;
         break;
      case PARSE_COMMENT:
         /* OMUPD DW 29/01/92 Just continue the while() loop */
         continue;
         break;
      case do_chain:
         fprintf(out,"Now Processing Degree of Freedom command \
Chain Closure:\n\n");
         p = (struct dof *) alloc(sizeof(struct dof));
         p->atoms = null;
         p->desc = (struct chain_closure_d *)
                   alloc(sizeof(struct chain_closure_d));
         p->dof_type = chain_closure_t;
         break;
      case do_forward:
         fprintf(out,"Now Processing Degree of Freedom command \
Forward (Backbone):\n\n");
         p = (struct dof *) alloc(sizeof(struct dof));
         p->atoms = null;
         p->desc = (struct backbone_d *) alloc(sizeof(struct backbone_d));
         p->dof_type = backbone_t;
         break;
      case do_reverse:
         fprintf(out,"Now Processing Degree of Freedom command \
Reverse (Backbone):\n\n");
         p = (struct dof *) alloc(sizeof(struct dof));
         p->atoms = null;
         p->desc = (struct backbone_d *) alloc(sizeof(struct backbone_d));
         p->dof_type = backbone_t;
         break;
      case do_side:
         fprintf(out,"Now Processing Degree of Freedom command \
Sidechain:\n\n");
         p = (struct dof *) alloc(sizeof(struct dof));
         p->atoms = null;
         p->desc = (struct sidechain_d *) alloc(sizeof(struct sidechain_d));
         p->dof_type = sidechain_t;
         break;
      case do_write:
         fprintf(out,"Now Processing Degree of Freedom command \
Write Coordinates:\n\n");
         if((write_fp = fopen(g_strparam[0],"w"))==NULL)
         {
            fprintf(out,"Unable to open output conformations file: %s\n",
                   g_strparam[0]);
            die();
         }
         p = (struct dof *) alloc(sizeof(struct dof));
         p->atoms = null;
         p->desc = (struct write_coordinates_d *)
                   alloc(sizeof(struct write_coordinates_d));
         p->dof_type = write_coordinates_t;
         break;
      case do_status:
         fprintf(out,"Now Processing Degree of Freedom command \
Write Status:\n\n");
         if((status_fp = fopen(g_strparam[0],"w"))==NULL)
         {
          fprintf(out,"Fatal: Unable to open status file: %s\n",g_strparam[0]);
             die();
         }
         p = (struct dof *) alloc(sizeof(struct dof));
         p->atoms = null;
         p->desc = (struct status_d *) alloc(sizeof(struct status_d));
         p->dof_type = status_t;
         break;
      case do_loops:
         fprintf(out,"Now Processing Degree of Freedom command \
Read Database Loops:\n\n");
         if((loops_fp = fopen(g_strparam[0],"r"))==NULL)
         {
             fprintf(out,"Fatal: Unable to open database loops file: %s\n",g_strparam[0]);
             die();
         }
         p = (struct dof *) alloc(sizeof(struct dof));
         p->atoms = null;
         p->desc = (struct rbest_d *) alloc(sizeof(struct rbest_d));
         p->dof_type = rbest_t;
         break;
 
#ifdef EVALUATE
      case do_evaluate:
         fprintf(out,"Now Processing Degree of Freedom command \
Evaluate:\n\n");
         p = (struct dof *) alloc(sizeof(struct dof));
         p->atoms = null;
         p->desc = (struct evaluate_d *) alloc(sizeof(struct evaluate_d));
         p->dof_type = evaluate_t;
         break;
#endif
 
      case do_end:
         fprintf(out,"\nEnd of commands. Conformational search will \
begin.\n\n");
         goto breakout;  /* If we don't do this, we'll add to our
                            linked list incorrectly */
         break;
 
      default:
         fprintf(out,"Command cannot follow CGEN: %s.\n",comline);
         die();
         break;
      }
 
      /* Insert this degree of freedom into the linked list. */
      p->next = p->prev = null;
      if (dof_head == null)
         dof_head = dof_tail = p;
      else
      {
         dof_tail->next = p;
         p->prev = dof_tail;
         dof_tail = p;
      }
      p->ok = true;
      ok = true;
 
      /* Now we setup other stuff for this degree of freedom
         calling appropriate routines to fill in the details
      */
      switch (key)
      {
      case do_chain:
         {  struct chain_closure_d *dp;
 
            dp = (struct chain_closure_d *)p->desc; /* ACRM added cast */
 
 
            setchncls(&dp->maxdt,&dp->maxg,&dp->maxevdw,
                                &dp->cistrans,&dp->startres,&ok);
            startres = dp->startres;
            lastres = startres+2;
            if (ok) chkchnclsr(startres,lastres,&ok);
            if (ok) finchncls(p,dp);
            p->ok = ok;
         }
         break;
      case do_forward:
         {  struct backbone_d *dp;
 
            dp = (struct backbone_d *)p->desc; /* ACRM added cast */
 
            dp->forward  = f77_true;
 
            sscanf(g_strparam[0],"%f",&dp->maxevdw);
 
            /* For the moment this is always true */
            dp->cistrans = f77_true;
 
            /* Used for N and C termini
               0.0 = use same as main map
            */
            dp->grid     = 0.0;
 
            /* These override calculation of the closing distance.
               if(qdelta) the closing distance is added to the
               calc'd value; otherwise it overrides
            */
            dp->qdelta = f77_false;
            dp->closing_distance = 0.0;
 
            /* Always true for now. If(false) symmetry won't be
               used for terminal residues
            */
            dp->tersym = f77_true;
 
            /* Max bond angle variation default is 5degrees
               PARSE_BACKBONE seems to set it to -1
            */
            dp->maxdt = -1.0;
 
            padspace(g_strparam[1],4);
            padspace(g_strparam[2],4);
            padspace(g_strparam[3],4);
            padspace(g_strparam[4],4);
            padspace(g_strparam[5],4);
 
            startres = get_resnum(g_strparam[1],g_strparam[2]);
            lastres  = get_resnum(g_strparam[1],g_strparam[3]);
            dp->closing_atom =
               get_atnum(g_strparam[1],g_strparam[4],g_strparam[5]);
 
            chkbbone(startres,lastres,&ok);
            propbbones(&p,dp,startres,lastres);
            p->ok = ok;
         }
         break;
 
      case do_reverse:
         {  struct backbone_d *dp;
 
            dp = (struct backbone_d *)p->desc; /* ACRM added cast */
 
            dp->forward  = f77_false;
 
            sscanf(g_strparam[0],"%f",&dp->maxevdw);
 
            /* For the moment this is always true */
            dp->cistrans = f77_true;
 
            /* Used for N and C termini
               0.0 = use same as main map
            */
            dp->grid     = 0.0;
 
            /* These override calculation of the closing distance.
               if(qdelta) the closing distance is added to the
               calc'd value; otherwise it overrides
            */
            dp->qdelta = f77_false;
            dp->closing_distance = 0.0;
 
            /* Always true for now. If(false) symmetry won't be
               used for terminal residues
            */
            dp->tersym = f77_true;
 
            /* Max bond angle variation default is 5degrees
               PARSE_BACKBONE seems to set it to -1
            */
            dp->maxdt = -1.0;
 
            padspace(g_strparam[1],4);
            padspace(g_strparam[2],4);
            padspace(g_strparam[3],4);
            padspace(g_strparam[4],4);
            padspace(g_strparam[5],4);
 
            startres = get_resnum(g_strparam[1],g_strparam[2]);
            lastres  = get_resnum(g_strparam[1],g_strparam[3]);
            dp->closing_atom =
               get_atnum(g_strparam[1],g_strparam[4],g_strparam[5]);
 
            chkbbone(startres,lastres,&ok);
            propbbones(&p,dp,startres,lastres);
            p->ok = ok;
         }
         break;
      case do_side:
         {  struct sidechain_d *dp;
            int *res;
            float *sgrid,*maxevdw;
            logical *vavoid,*clump_symmetry;
            int evalopt;
 
            dp = (struct sidechain_d *)p->desc; /* ACRM added cast */
            res = (int *) alloc(sizeof(int)*(values.nres+1));
            sgrid = (float *) alloc(sizeof(float)*(values.nres+1));
            maxevdw = (float *) alloc(sizeof(float)*(values.nres+1));
            vavoid = (logical *) alloc(sizeof(logical)*(values.nres+1));
            clump_symmetry = (logical *) alloc(sizeof(logical)*(values.nres+1));
 
            setsidchn(startres,lastres,dp,
                            res,sgrid,maxevdw,vavoid,clump_symmetry,&ok);
 
            finsidchn(p,dp,res,sgrid,maxevdw,vavoid,clump_symmetry);
            free(res);
            free(sgrid);
            free(maxevdw);
            free(vavoid);
            free(clump_symmetry);
            startres = 0;
            lastres = 0;
            p->ok = ok;
         }
         break;
      case do_write:
         {  struct write_coordinates_d *dp;
 
            dp = (struct write_coordinates_d *)p->desc; /* ACRM added cast */
            dp->maxconf = largint;
            dp->use_comp = f77_true;
            p->ok = f77_true;
         }
         break;
      case do_status:
         /* ACRM added cast */
         setup_status((struct status_d *)p->desc,&ok);
         break;
      case do_loops:
         {  struct rbest_d *dp;
 
            dp = (struct rbest_d *)p->desc; /* ACRM added cast */
            dp->nbest = 1000; /* Doesn't matter since all energies are
                                 the same so all will be read */
            dp->maxevdw = largnum;
            gtrbstatms(p,dp,&ok);
            p->ok = f77_true;
         }
         break;
 
#ifdef EVALUATE
      case do_evaluate:
         {
            struct evaluate_d *dp;
 
            dp = p->desc;
            dp->minist = alloc(8+sizeof(char));
            strcpy(dp->minist,"ENERGY");
            dp->ministl = strlen(dp->minist);
            dp->eval_code = eval_code_energy;
            p->ok = f77_true;
         }
         break;
#endif
 
      default:
         break; /* We should have intercepted errors already */
      }
      if (!ok) parse_ok = false;
   }
 
/* Jump out from while() loop if the END command has been reached */
breakout:
 
   /* Sets the grid limits */
   extnlims();
 
/***   We don't handle titles
//   rdtitl(cg.ctitle,&cg.nctitl,&five);
//   fprintf(out,"\nTitle for written conformations:\n\n");
//   prtitl(cg.ctitle,&cg.nctitl,&six);
***/
 
   setrestitr();
   prtcgncmnd();
   if (!parse_ok)
      fprintf(out,"Due to errors in parsing, CGEN will not be run\n");
   return parse_ok;
}
 

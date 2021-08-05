#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
#include "cgparse.h" 

void prtcgncmnd(
void
)
{
   struct dof *p;
   int count,iseg,ires;
 
   fprintf(out,"\n\nResults of parsing the CGEN commands:\n\n");
   prtglobopt();
   for  (p=dof_head,count=1; p; p=p->next,count++)
   {
      fprintf(out,"\nDegree of freedom %d is %s.  Restart iteration is %d\n",
                 count,move_name[p->dof_type],p->restart_iter);
      if (!p->ok)
         fprintf(out,"  Parsing ERRORS were found for this degree of freedom.\n");
      else
      {
         switch (p->dof_type)
         {
         case chain_closure_t:
            {
               struct chain_closure_d *dp;
 
               dp = (struct chain_closure_d *)p->desc; /* ACRM added cast */
               iseg = getseg(&dp->startres,pstruct.segndx,&values.nsegs) - 1;
               ires = dp->startres - 1;
               fprintf(out,"Segment: %.4s  First Residue ID: %.4s  \
Residue names: %.4s %.4s %.4s\n", 
                       pstruct.segid[iseg],pstruct.resid[ires],
                       pstruct.resnme[ires],pstruct.resnme[ires+1],
                       pstruct.resnme[ires+2]);
               fprintf(out,"MAXDT = %g  MAXG = %g\n",
                      dp->maxdt/dtorad,dp->maxg);
               fprintf(out,"Cis-trans isomerization will%s be considered.\n",
                      dp->cistrans ? "" : " not");
               print_atoms(p->atoms);
            }
            break;
         case backbone_t:
            {
               struct backbone_d *dp;
 
               char buf[21];
               dp = (struct backbone_d *)p->desc; /* ACRM added cast */
               iseg = getseg(&dp->resno,pstruct.segndx,&values.nsegs) - 1;
               ires = dp->resno - 1;
               fprintf(out,"Segment: %.4s  Residue ID: %.4s  Residue name: %.4s\n",

                       pstruct.segid[iseg],pstruct.resid[ires],
                       pstruct.resnme[ires]);
               fprintf(out,"MAXEVDW = %g\n",dp->maxevdw);
               fprintf(out,"Cis-trans isomerization will%s be considered.\n",
                      dp->cistrans ? "" : " not");
               fprintf(out,"Construction will be in the %s direction.\n",
                      dp->forward ? "N to C" : "C to N");
               if (dp->closing_atom == 0)
                   fprintf(out,"There is no closing atom.\n");
               else
               {
                  buf[0] = '\0';
                  print_atom1(buf,dp->closing_atom);
                  fprintf(out,"Closing atom is %s and \ncurrent closing distance \
is %s%g.\n", buf, dp->qdelta ? "Delta " : "", dp->closing_distance);
                  if (dp->maxdt == -1.0)
                  {
                     fprintf(out,"MAXDT will be determined from the degrees \
of freedom.\n");
                  }
                  else
                  {
                     fprintf(out,"MAXDT = %g\n",dp->maxdt/dtorad);
                  }
               }
               if (dp->nter || dp->cter)
               {
                  fprintf(out,"This residue is at the %s terminus\n",
                         dp->nter ? "N" : "C");
                  fprintf(out,"The backbone grid is currently %g with %ssymmetry.\n",
                         dp->grid,dp->tersym ? "" : "no ");
               }
               print_atoms(p->atoms);
            }
            break;
         case sidechain_t:
            {
               struct sidechain_d *dp;
               dp = (struct sidechain_d *)p->desc; /* ACRM added cast */
               fprintf(out,"Sidechain construction option is %s\n",
                      st_sideopt[dp->sideopt-1]);
               fprintf(out,"Evaluation option is %s\n",st_evalopt[dp->evalopt-1]);
               fprintf(out,"NCOMB = %d  MAXSIDEITER = %d\n",
                      dp->ncomb,dp->maxsideiter);
               prtsidchn(dp);
            }
            break;
         case write_coordinates_t:
            {
               struct write_coordinates_d *dp;
               dp = (struct write_coordinates_d *)p->desc; /* ACRM added cast */
               fprintf(out,"CUNIT = %d  MAXCONF = %d\n",dp->cunit,dp->maxconf);
               if (dp->use_comp)
                  fprintf(out,"Comparison coordinates will be used for \
reference.\n");
            }
            break;
         case rbest_t:
            {
               struct rbest_d *dp;
               dp = (struct rbest_d *)p->desc; /* ACRM added cast */
               fprintf(out,"UNIT = %d  NBEST = %d  MAXEVDW = %g\n",
                      dp->unit,dp->nbest,dp->maxevdw);
               print_atoms(p->atoms);
            }
            break;
 
#ifdef EVALUATE
         case evaluate_t:
            {
               struct evaluate_d *dp;
               dp = p->desc;
               dp->minist[dp->ministl] = '\0';
               switch (dp->eval_code)
               {
               case eval_code_rms:
                  fprintf(out,"Evaluating RMS deviation from starting \
coordinates\n");
                  break;
               case eval_code_energy:
                  fprintf(out,"Evaluating energy using minimization command:\n%s\n",
                         dp->minist);
                  break;
               case eval_code_user:
                  fprintf(out,"Evaluation by calling the user subroutine, \
user_cgeval.\n");
                  break;
               }
            }
            break;
#endif
 
         case status_t:
            {
               struct status_d *dp;
               int i;
               dp = (struct status_d *)p->desc; /* ACRM added cast */
               if (dp->setprn)
                  fprintf(out,"Process name will be set to status.\n");
               if (dp->fileunit >=0)
               {
                  fprintf(out,"Status information will be written to unit %d \
every %d calls.\n", dp->fileunit,dp->filefreq);
               }
               if (dp->nflush > 0)
               {
                  fprintf(out,"The following output units will be flushed:\n");
                  for (i=0; i<dp->nflush; i++)
                       fprintf(out,"Unit %d every %d calls",
                              dp->flushunit[i],dp->flushfreq[i]);
                  fprintf(out,"and every hour.\n");
               }
            }
            break;
         }
      }
   }
}
 

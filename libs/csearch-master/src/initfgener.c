#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
void initfgener(
void
)
{
struct dof *p;
char title[80];
int l,i;
 
   rdpepmaps();
   rdprocons(cg.proconsu,&cg.procons,&cg.proconsphi,
             &cg.eprocons,&cg.nprocons,cg.eringpro);
   getgrdspc();

   fill_grid();

   setup_nbond();
   setup_hbond();
   if (cg.save_coor)
   {
/*  Following commented out until HBOND can be stopped from generating
    hydrogen positions.
      setup_consp(&cg.consp,&cg.nconsp,(struct dof *)null);
      cg.savex = (float *) alloc(cg.nconsp * sizeof(float));
      cg.savey = (float *) alloc(cg.nconsp * sizeof(float));
      cg.savez = (float *) alloc(cg.nconsp * sizeof(float));
      for(i=0; i<cg.nconsp; i++)
      {
         i1 = cg.consp[i] - 1;
         cg.savex[i] = coords.xcart[i1];
         cg.savey[i] = coords.ycart[i1];
         cg.savez[i] = coords.zcart[i1];
      }
*/
      cg.savex = (float *) alloc(values.natoms * sizeof(float));
      cg.savey = (float *) alloc(values.natoms * sizeof(float));
      cg.savez = (float *) alloc(values.natoms * sizeof(float));
      CopyFArray(coords.xcart,cg.savex,values.natoms);
      CopyFArray(coords.ycart,cg.savey,values.natoms);
      CopyFArray(coords.zcart,cg.savez,values.natoms);
   }
   for (p=dof_head; p; p=p->next)
   {
      switch (p->dof_type)
      {
      case write_coordinates_t:
         {
            struct write_coordinates_d *dp;
 
            dp = (struct write_coordinates_d *)p->desc; /* ACRM added cast */
            setup_consp(&dp->consp,&dp->nconsp,p->prev);
/*
//            cgwrtinit(&dp->cunit,cg.ctitle,&cg.nctitl,dp->consp,
//                          &dp->nconsp);
*/
            /* C version */
            cgwrtinit(dp->consp,dp->nconsp);
            strcpy(title,"Starting coordinates before CGEN");
            filspc(title,&eighty,(l=strlen(title),&l));
 
/* We don't allow comparison sets
//            if (dp->use_comp)
//            {
//              xbyind(sizeof coords.xcart[0],coords.xcart,coords.xwork,
//                             dp->consp,dp->nconsp);
//               xbyind(sizeof coords.ycart[0],coords.ycart,coords.ywork,
//                             dp->consp,dp->nconsp);
//               xbyind(sizeof coords.zcart[0],coords.zcart,coords.zwork,
//                             dp->consp,dp->nconsp);
//            }
*/
/*
//            cgwrtnxt(&dp->cunit,dp->consp,&dp->nconsp,title,&one,&fzero);
*/
            /* C version */
            cgwrtnxt(dp->consp,dp->nconsp,0.0);
            if (dp->use_comp)
            {
               xbyind(sizeof coords.xcart[0],(char *)coords.xcart,(char *)coords.xwork,
                             dp->consp,dp->nconsp);
               xbyind(sizeof coords.ycart[0],(char *)coords.ycart,(char *)coords.ywork,
                             dp->consp,dp->nconsp);
               xbyind(sizeof coords.zcart[0],(char *)coords.zcart,(char *)coords.zwork,
                             dp->consp,dp->nconsp);
            }
         }
         break;
      case chain_closure_t:
/* ACRM added cast */
         gtclschngom((struct chain_closure_d *)p->desc);
         break;
#ifdef EVALUATE
      case evaluate_t:
         {
            struct evaluate_d *dp;
            dp = p->desc;
            dp->evalcount = 0;
            setup_consp(&dp->consp,&dp->nconsp,p->prev);
            setup_imove(dp);
            dp->workst = alloc(dp->ministl);
/*
*   The following code is commented out until HBOND stops generating
*   hydrogen positions.
            dp->savedx = (float *) alloc(dp->nconsp * sizeof(float));
            dp->savedy = (float *) alloc(dp->nconsp * sizeof(float));
            dp->savedz = (float *) alloc(dp->nconsp * sizeof(float));
*/
            dp->savedx = (float *) alloc(values.natoms * sizeof(float));
            dp->savedy = (float *) alloc(values.natoms * sizeof(float));
            dp->savedz = (float *) alloc(values.natoms * sizeof(float));
            dp->initx = (float *) alloc(dp->nconsp * sizeof(float));
            dp->inity = (float *) alloc(dp->nconsp * sizeof(float));
            dp->initz = (float *) alloc(dp->nconsp * sizeof(float));
            for (i=0; i < dp->nconsp; i++)
            {
               dp->initx[i] = coords.xcart[dp->consp[i]-1];
               dp->inity[i] = coords.ycart[dp->consp[i]-1];
               dp->initz[i] = coords.zcart[dp->consp[i]-1];
            }
         }
         break;
#endif   /* End of #ifdef EVALUATE */
      case status_t:
         {
            struct status_d *dp;
 
            dp = (struct status_d *)p->desc; /* ACRM added cast */
            dp->callcount = 0;
	    /* OMUPD bnj 20/11/91 correct arg list to funct time()... */
            dp->file_time = time(NULL);
            for (i=0; i<dp->nflush; i++) dp->flush_time[i] = time(NULL);
            dp->file_count = 0;
         }
         break;
      }
   }
/*
*   The completion of the closing atom data structure must be done after
*   the chain closures have been completed.
*/
   for (p=dof_head; p; p=p->next)
   {
      switch (p->dof_type)
      {
         case backbone_t:
/* ACRM added cast */
            finclsatm((struct backbone_d *)p->desc);
/* ACRM added cast */
            chkftres((struct backbone_d *)p->desc);
            break;
      }
   }
   walk_atoms(int i1;
              i1 = atp->atomno - 1;
              coords.xcart[i1] = coords.ycart[i1] = coords.zcart[i1] = anum;);
   dummy_sidehits = (int *) alloc(values.natoms * sizeof(int));
}




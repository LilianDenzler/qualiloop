#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Completes the construction of the sidechain degree of freedom
   which basically entails filling in all the fields given the list of
   residues stored in res. res is zero terminated.
*/
 
void finsidchn(
struct dof *dofp,
struct sidechain_d *desc,
int *res,
float *sgrid,
float *maxevdw,
logical *vavoid,
logical *clump_symmetry
)
{
   int ires,nr,istres,istclump,istatm,count;
   struct sideres *srp,**srpp;
   struct clump *clp,**clpp;
   struct atom *atp,**atpp,**dof_atpp;
   int i,j,k,parmcode,icphi;
/*
   double sidchnengy();
   double sidchnrms();
*/
   int *ip;
   float *xp,*yp,*zp;
 
   switch (desc->evalopt)
   {
   case eo_energy:
      /*desc->eval_func = &sidchnengy;*/
      desc->eval_func = sidchnengy;
      break;
   case eo_rms:
      /*desc->eval_func = &sidchnrms;*/
      desc->eval_func = sidchnrms;
      break;
   default:
      fprintf(out,"\nError in finsidchn -- Bad evalopt = %d\n",
              desc->evalopt);
      die();
   }
   desc->energy = 0.0;
   for (nr=0; res[nr]; nr++);
   count = 0;
   desc->residues = (struct sideres **) alloc((nr+1)*sizeof(struct sideres *));
   for (srpp=desc->residues; (ires = *res); res++, srpp++)
   {
/* ACRM Changed to use C version of srchwd()
//      istres = srchwd(sidetop.sc_resname,&sidetop.nscres,&pstruct.resnme[ires-1]);
*/
      istres = srchwd(sidetop.sc_resname,sidetop.nscres,pstruct.resnme[ires-1]);
      if (istres == 0)
      {
         fprintf(out,"\nError in finsidchn -- Unable to find %.4s in \
sidechain topology\n", pstruct.resnme[ires-1]);
         die();
      }
      srp = *srpp = (struct sideres *) alloc(sizeof(struct sideres));
      srp->sidehits = (int *) alloc(values.natoms * sizeof(int));
      srp->special = sidetop.sc_special[istres-1];
      srp->nclump = sidetop.sc_clump_part[istres] -
                    sidetop.sc_clump_part[istres-1];
      srp->clumps = (struct clump **)
                    alloc((srp->nclump+1) * sizeof(struct clump *));
      srp->resno = ires;
      srp->sgrid = *sgrid++;
      srp->maxevdw = *maxevdw++;
      srp->vavoid = *vavoid++;
      srp->clump_symmetry = *clump_symmetry++;
      if (srp->vavoid && srp->maxevdw < 0.0)
      {
         fprintf(out,"\nWarning from FINISH_SIDECHAIN for sidechain %.4s %.4s\n",
                pstruct.resnme[srp->resno-1],pstruct.resid[srp->resno-1]);
         fprintf(out,"Maxevdw must not be negative when van der Waals avoidance \
is in effect.\n");
         fprintf(out,"Maxevdw will be set to zero.\n");
         srp->maxevdw = 0.0;
      }
      if (srp->sgrid > 0.0) srp->sgrid *= dtorad;
      clpp = srp->clumps;
      for (istclump=sidetop.sc_clump_part[istres-1]+1;
          istclump<=sidetop.sc_clump_part[istres];
          istclump++,clpp++)
      {
         clp = *clpp = (struct clump *) alloc(sizeof(struct clump));
         clp->symmetry = srp->clump_symmetry ?
                         sidetop.sc_symmetry[istclump-1] : 1;
         clp->natom = sidetop.sc_atom_part[istclump] -
                      sidetop.sc_atom_part[istclump-1];
         clp->best_phi = (float *) alloc(desc->ncomb * sizeof(float));
         atpp = clp->atoms = (struct atom **)
                             alloc(sizeof(struct atom *) * (clp->natom+1));
         for (istatm=sidetop.sc_atom_part[istclump-1]+1;
              istatm<=sidetop.sc_atom_part[istclump];
              istatm++,atpp++)
         {
            atp = *atpp = (struct atom *) alloc(sizeof(struct atom));
            count++;
            assign_name(pstruct.atmnme,pstruct.lstatm,&ires,&atp->ante[0],
                        sidetop.sc_ante3_bld[istatm-1]);
            assign_name(pstruct.atmnme,pstruct.lstatm,&ires,&atp->ante[1],
                        sidetop.sc_ante2_bld[istatm-1]);
            assign_name(pstruct.atmnme,pstruct.lstatm,&ires,&atp->ante[2],
                        sidetop.sc_ante1_bld[istatm-1]);
            assign_name(pstruct.atmnme,pstruct.lstatm,&ires,&atp->atomno,
                        sidetop.sc_atom_bld[istatm-1]);
            atp->rcontact = null;
            atp->bond = sidetop.sc_bond_bld[istatm-1];
            if (!atp->bond)
                 atp->bond = getparbond2(atp->ante[2],atp->atomno);
            atp->angle = sidetop.sc_angle_bld[istatm-1];
            if (!atp->angle)
                 atp->angle = gtprangl2(atp->ante[1],atp->ante[2],
                                           atp->atomno);
            atp->torsion = sidetop.sc_tors_bld[istatm-1];
            atp->offset = sidetop.sc_offset[istatm-1];
            atp->cons_code = sidetop.sc_code_bld[istatm-1];
            if (atp->cons_code == st_free)
            {
               j = pstruct.atcode[atp->ante[1]-1];
               k = pstruct.atcode[atp->ante[2]-1];
               parmcode = (int)(*(cg.parm_no + 100*(k-1) + j-1)); /* ACRM 20.02.06 FIXME */
               icphi = nindx(&parmcode,engpar.torkey,&values.nptpar);
               if (icphi == 0)
               {
                  fprintf(out,"\nError in finsidchn -- No \
parameters for torsion between %.4s and %.4s:\n",
                             pstruct.atmnme[atp->ante[1]-1],
                             pstruct.atmnme[atp->ante[2]-1]);
                  die();
               }
               atp->torsion = (-pi-engpar.torphs[icphi-1]) / 
                                engpar.tormlt[icphi-1];
               atp->torsion_period = 2 * pi / engpar.tormlt[icphi-1];
            }
            if (strncmp(sidetop.sc_atom_bld[istatm-1],
                        sidetop.sc_free_atom[istclump-1],4) == 0)
            clp->free_atom = atp;
         }
         *atpp = null;
         if (srp->vavoid) setavdclp(srp,clp);
      }
      *clpp = null;
      srp->beste = (float *) alloc(desc->ncomb * sizeof(float));
      srp->nconf = 0;
      srp->natom = 0;
      walk_over_sideres(srp,srp->natom++;);
      ip = srp->atomp = (int *) alloc((srp->natom+1)*sizeof(int));
      walk_over_sideres(srp,*ip++ = atp->atomno;);
      *ip = 0;
      srp->bestxpp = (float **) alloc(desc->ncomb * sizeof(float *));
      srp->bestypp = (float **) alloc(desc->ncomb * sizeof(float *));
      srp->bestzpp = (float **) alloc(desc->ncomb * sizeof(float *));
      for (i=0; i<desc->ncomb; i++)
      {
         srp->bestxpp[i] = (float *) alloc(srp->natom * sizeof(float));
         srp->bestypp[i] = (float *) alloc(srp->natom * sizeof(float));
         srp->bestzpp[i] = (float *) alloc(srp->natom * sizeof(float));
      }
      xp = srp->refx = (float *) alloc(srp->natom * sizeof(float));
      yp = srp->refy = (float *) alloc(srp->natom * sizeof(float));
      zp = srp->refz = (float *) alloc(srp->natom * sizeof(float));
      for (ip=srp->atomp; (i = *ip); ip++)
      {
         *xp++ = coords.xcart[--i];
         *yp++ = coords.ycart[i];
         *zp++ = coords.zcart[i];
      }
      selsidetor(srp,&srp->torsions,pstruct.attor1,pstruct.attor2,
                 pstruct.attor3,pstruct.attor4,getpar.itorp,values.nptors);
      selsidetor(srp,&srp->impropers,pstruct.atimp1,pstruct.atimp2,
                 pstruct.atimp3,pstruct.atimp4,getpar.iimpp,values.nitors);
   }
   *srpp = null;
   dofp->atoms = dof_atpp =
      (struct atom **) alloc((count+1)*sizeof(struct atom *));
   desc->nresidues = 0;
   for (srpp=desc->residues; *srpp; srpp++)
   {
      desc->nresidues++;
      for (clpp=(*srpp)->clumps; *clpp; clpp++)
         for (atpp=(*clpp)->atoms; *atpp; atpp++)
            *dof_atpp++ = *atpp;
   }
   *dof_atpp = null;
}
 

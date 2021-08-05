#include "ProtoTypes.h"
#include "CongenProto.h"
 
#define MAXCLS 65
 
void chnclsr_f(
struct dof *dofp,
struct chain_closure_d *desc
)
{
   int ires,i,ind,i1;
   float newx[MAXCLS][6],newy[MAXCLS][6],newz[MAXCLS][6],
         clsomega[MAXCLS][6],pepomega[MAXCLS][3],newangle[MAXCLS][3][3],psi;
   logical ok,badphi;
   short short1,short2,short3,short4;
   int iomega;
   double sum,dx,dy,dz;
   struct atom *atp;
   int iter1,ncls;
 
   iter1 = 0;
   ncls = 0;
   do
   {
      if (ncls > 0)
      {
         /* ACRM This is a major hack! FIXME!!! */
         CopyFArray((float *)(newangle[ncls-1]),(float *)(newangle[ncls]),9);

         /* These are OK */
         CopyFArray(clsomega[ncls-1],clsomega[ncls],6);
         CopyFArray(pepomega[ncls-1],pepomega[ncls],3);
      }
      /* COMUPD drm 20/11/91 change name of routine to add a */
      clschna(coords.xcart, coords.ycart, coords.zcart,
             desc->atmind, desc->clsbond, desc->clsangle,
             newx[ncls], newy[ncls], newz[ncls],
             clsomega[ncls], &iter1, desc->clscistrans,
             pepomega[ncls], &desc->maxdt, newangle[ncls], &desc->maxg);
 
      if (iter1 > 0)
      {
         if (++ncls >= MAXCLS)
            fprintf(out,"Error in chnclsr_f -- MAXCLS (%d) exceeded.\n",
                    ncls);
      }
   }
 
   while (iter1 > 0);
 
   if (dbg.cgen) fprintf(out,"Chain_closure_f: %d chain closures found.\n",
                         ncls);
   if (dofp->restart_iter > 1)
   {
      dofp->iter = dofp->restart_iter;
      dofp->restart_iter = 0;
   }
   else
      dofp->iter = 1;
 
   for (; dofp->iter<=ncls; dofp->iter++)
   {
      if (dbg.cgen) fprintf(out,"Chain_closure_f: Iteration %d\n",dofp->iter);
      iter1 = dofp->iter-1;
      ok = true;
      for (ires=desc->startres; ires<=desc->startres+2; ires++)
      {
         ind = (ires - desc->startres + 1)*2 - 1;
         badphi = strncmp("PRO ",pstruct.resnme[ires-1],4) == 0 &&
                  (clsomega[iter1][ind-1] < cg.proconsphi[0] ||
                  clsomega[iter1][ind-1] > cg.proconsphi[cg.nprocons-1]);
         if (badphi)
         {
            if (dbg.cgen > 1)
               fprintf(out,"Proline phi angle %g is out of range.\n",
                          clsomega[iter1][ind-1]/dtorad);
            ok = false;
         }
      }
      if (ok)
      {
         atp = desc->first_n;
         phia((short1 = atp->ante[0],&short1),
              (short2 = atp->ante[1],&short2),
              (short3 = atp->ante[2],&short3),
              (short4 = atp->atomno,&short4),
              &psi,
              &one,
              coords.xcart,coords.ycart,coords.zcart);
         psi *= dtorad;
         delatmfgrd(atp->atomno);
         atp->bond = desc->clsbond[0][1];
         atp->angle = newangle[iter1][0][0];
         constatm(atp,psi);
         adatmtgrd(atp->atomno);
         for (i=0,iomega=0; i<=2; i++)
         {
            if (i < 2)
            {
               desc->at_ca[i]->angle = newangle[iter1][i][1];
               desc->at_c[i]->angle = newangle[iter1][i][2];
               desc->at_n[i]->angle = newangle[iter1][i+1][0];
            }
            constatm(desc->at_ca[i],pepomega[iter1][i]);
            constatm(desc->at_h[i],pepomega[iter1][i]-pi);
            constatm(desc->at_c[i],clsomega[iter1][iomega]);
            if (strncmp(pstruct.resnme[desc->startres+i-1],"PRO ",4) == 0)
               mkprolrng(coords.xcart,coords.ycart,coords.zcart,
                                 &desc->at_cb[i]->ante[0],
                                 &desc->at_cb[i]->ante[1],
                                 &desc->at_cb[i]->ante[2],
                                 &desc->at_cb[i]->atomno,
                                 &desc->at_cg[i]->atomno,
                                 &desc->at_cd[i]->atomno,
                                 &clsomega[iter1][iomega++],
                                 cg.procons,cg.proconsphi,&cg.nprocons);
            else
               constatm(desc->at_cb[i],clsomega[iter1][iomega++]-rad120);
 
            constatm(desc->at_n[i],clsomega[iter1][iomega]);
            constatm(desc->at_o[i],clsomega[iter1][iomega++]-pi);
         }
         for (sum=0.0,i=1; i<=6; i++)
         {
            i1 = desc->atmind[i][0] - 1;
            dx = coords.xcart[i1] - newx[iter1][i-1];
            dy = coords.ycart[i1] - newy[iter1][i-1];
            dz = coords.zcart[i1] - newz[iter1][i-1];
            sum += dx*dx + dy*dy + dz*dz;
         }
         if (sqrt(sum/6.0) > 1.0e-4)
            fprintf(out,"Warning from chnclsr_f -- CARTX2 and CLSCHN \
RMS deviation is %.3g\n",
                       sqrt(sum/6.0));
         if (chkcntcts(dofp->atoms,&desc->maxevdw,dummy_sidehits))
         {
            dispatch(dofp->next);
            delatmsfgrd(dofp->atoms);
         }
      }
   }
}
 

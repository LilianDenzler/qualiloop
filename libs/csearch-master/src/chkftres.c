#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/*
*   Checks if this backbone degree of freedom is for a terminal
*   residue, and if so, checks if the terminal grid is defined.
*   If not, then it is set to the minimum grid separation found
*   in the peptide maps.
*
*   Once the grid size is in hand, a grid for a terminal residue is
*   constructed if necessary.
*/
 
void chkftres(
struct backbone_d *desc
)
{
   int interval,size;
   struct pepmap *map;
   char *res;
   int nmap,i,j,iseg;
   float mingrid,d1,d2,d3;
   float omega_limit,omega_inc;
   float phi_limit,phi_inc;
   float psi_limit,psi_inc;
   float omega,phi,psi;
 
   iseg = getseg(&desc->resno,pstruct.segndx,&values.nsegs);
   desc->ntermap = 0;
   desc->termap = null;
   desc->grid *= dtorad;
   if (desc->nter || desc->cter)
   {
      if (desc->nter && strncmp(pstruct.resnme[desc->resno-1],"PRO ",4) == 0 &&
          desc->tersym)
      {
         desc->tersym = false;
         fprintf(out,"Terminal symmetry being turned off for residue \
%.4s %.4s %.4s\n", pstruct.segid[iseg-1],
                pstruct.resnme[desc->resno-1],
                pstruct.resid[desc->resno-1]);
      }
      if (desc->grid == 0.0)
      {
         /*res = &pstruct.resnme[desc->resno-1]; ACRM & may not be needed here */
/* COM BNJ UPDOML remove address operator... */

         res = pstruct.resnme[desc->resno-1];  
         if (strncmp(res,"PRO ",4) == 0)
         {
            map = cg.promap;
            nmap = cg.npromap;
         }
         else if (strncmp(res,"GLY ",4) == 0)
         {
            map = cg.glymap;
            nmap = cg.nglymap;
         }
         else
         {
            map = cg.aamap;
            nmap = cg.naamap;
         }
         mingrid = largnum;
         for (i=0; i<nmap; i++)
         {
            for (j=0; j<i; j++)
            {
               d1 = fabs((map+i)->omega - (map+j)->omega);
               d2 = fabs((map+i)->phi - (map+j)->phi);
               d3 = fabs((map+i)->psi - (map+j)->psi);
               if (d1 == 0.0) d1 = largnum;
               if (d2 == 0.0) d2 = largnum;
               if (d3 == 0.0) d3 = largnum;
               mingrid = minf3(d1,d2,minf2(d3,mingrid));
            }
         }
         desc->grid = mingrid;
         fprintf(out,"For residue %.4s %.4s %.4s, the terminal residue grid \
is set to %g degrees.\n",
                pstruct.segid[iseg-1], res,
                pstruct.resid[desc->resno-1],
                desc->grid/dtorad);
      }
      omega_inc = pi;
      omega_limit = pi/2.0;
      phi_limit = psi_limit = pi - desc->grid/2.0;
      phi_inc = psi_inc = desc->grid;
      interval = 2*pi / desc->grid;
      size = interval * interval * 2;
      if (desc->nter)
      {
         omega_limit = -pi;
         if (desc->forward)
         {
            size = interval * interval;
         }
         else
         {
            if (desc->tersym)
            {
               size = (interval+2) / 3;
               phi_limit = maxf2(-pi/3.0 - desc->grid/2.0, -pi);
            }
            else size = interval;
            size *= interval;
         }
      }
      else
      { /* desc->cter */
         if (desc->forward && desc->tersym)
         {
            size = interval * interval;
            psi_limit = maxf2(-desc->grid/2.0, -pi);
         }
      }
      desc->ntermap = size;
      desc->termap = (struct pepmap *) alloc(sizeof(struct pepmap) * size);
      map = desc->termap;
      for (omega = -pi; omega <= omega_limit; omega += omega_inc)
      {
         for (phi = -pi; phi <= phi_limit; phi += phi_inc)
         {
            for (psi = -pi; psi <= psi_limit; psi += psi_inc)
            {
               map->omega = omega;
               map->phi = phi;
               map->psi = psi;
               map->energy = 0.0;
               map->cisonly = fabs(omega/dtorad) < 90.0;
               if (map++ > desc->termap+size)
               {
                  fprintf(out,"Error in chkftres -- Map (%lu) of \
size (%d) overrun.\n",(unsigned long)desc->termap,size);
                  die();
               }
            }
         }
      }
   }
}
 

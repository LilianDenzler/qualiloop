#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   This procedure maps the grid points specified through the bounds and
*   grid size into the range of allowed torsion angles in rangep. The
*   mapping procedure currently employed check each grid point to see if
*   it is allowed. If not, it is mapped to the closest allowed boundary
*   The code is complicated because the torsion angle space is periodic.
*   We assume that the range list contains elements only in [-pi,pi].
*/
 
void mapgrdtorng(
float lbound,
float hbound,
float grid,
struct range_list *rangep,
float **philistp,
int *nphip
)
{
   int nphi,count;
   float *phip,phi,d1,d2,nextl,cphi,*phi2p,*newphip,f,clow,chigh,
         phi1,phi2,dphi;
   int irp,nr;
   logical found;
 
   /* ACRM added cast */
   nphi = chmceil((f=(float)((hbound-lbound)/grid),&f))+1;
   phip = *philistp = (float *) alloc(sizeof(float)*(nphi+1));
   nr = rangep->nrange;
   if (nr == 0)
   {
      if (dbg.cgen > 1)
         fprintf(out,"Map_grid_into_range found no points\n");
      *nphip = 0;
 
return;
   }
   irp = 0;
   clow = *(rangep->low);
   chigh = *(rangep->high);
   for (phi=lbound; phi<=hbound; phi += grid)
   {
      if (phi < -pi)
         cphi = phi + 2*pi;
      else if
         (phi > pi) cphi = phi - 2*pi;
      else
         cphi = phi;
      found = false;
      count = 0;
      while (!found)
      {
         if (clow <= cphi && cphi <= chigh)
         {
            *phip++ = cphi;
            found = true;
         }
         else
         {
            nextl = rangep->low[irp+1<nr ? irp+1 : 0];
            if (nextl <= chigh)
            {
               if (chigh < cphi && cphi <= pi)
               {
                  found = true;
                  d1 = cphi - chigh;
                  d2 = nextl + 2*pi - cphi;
                  *phip++ = d1 < d2 ? chigh : nextl;
               }
               else if (-pi <= cphi && cphi < nextl)
               {
                  found = true;
                  d1 = cphi - chigh + 2*pi;
                  d2 = nextl - cphi;
                  *phip++ = d1 < d2 ? chigh : nextl;
               }
            }
            else if (chigh < cphi && cphi < nextl)
            {
               found = true;
               d1 = cphi - chigh;
               d2 = nextl - cphi;
               *phip++ = d1 < d2 ? chigh : nextl;
            }
         }
         if (!found)
         {
            if (count++ > nr)
            {
               fprintf(out,"Error in mapgrdtorng -- \
Infinite loop detected.\n");
               die();
            }
            if (++irp >= nr) irp = 0;
            clow = rangep->low[irp];
            chigh = rangep->high[irp];
         }
      }
   }
   nphi = phip - *philistp;
   if (dbg.cgen > 1)
   {
      int i,j,l;
      char buf[81];
      fprintf(out,"Mapped grid points are \n");
      for (i=0; i<nphi; i += 10)
      {
         buf[0] = '\0';
         l = mini2(i+9,nphi-1);
         for (j=i; j<=l; j++)
         {
            sprintf(&buf[(j-i)*8],"%7.3f ",(*philistp)[j]);
         }
         fprintf(out,"%s\n",buf);
      }
   }
   for (newphip = phi2p = *philistp; phi2p < phip; phi2p++)
   {
      if (phi2p == *philistp)
         *newphip++ = *phi2p;
      else
      {
         phi1 = *phi2p;
         phi2 = *(phi2p-1);
         dphi = fabs(phi1-phi2);
         if (!(dphi < 3.0*range_eps || fabs(dphi - 2*pi) < 3.0*range_eps))
            *newphip++ = *phi2p;
      }
   }
   *nphip = newphip - *philistp;
/*
*   Take care of circular lists.
*/
   if (*nphip > 1)
   {
      phi1 = *(newphip-1);
      phi2 = **philistp;
      dphi = fabs(phi1-phi2);
      if (dphi < 3.0*range_eps || fabs(dphi - 2*pi) < 3.0*range_eps)
         (*nphip)--;
   }
   if (dbg.cgen > 1)
      fprintf(out,"A total of %d points were used.\n",*nphip);
}
 

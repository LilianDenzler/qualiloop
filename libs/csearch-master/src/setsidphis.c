#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Determines the torsion angles to be used for sampling the sidechain
*   grid. The array of torsion angles is returned in philistp, which is
*   allocated by this routine, and which be freed by the calling
*   routine. If van der Waals repulsion avoidance is in effect then the
*   sidehits array will return which sidechains influenced the grid
*   selection.
*
*   The philist array is allocated with enough space to accomodate one
*   extra entry for use by the never_fail option in the sidchnbst
*   routine.
*/
 
void setsidphis(
struct sideres *srp,
struct clump *clp,
float **philistp,
int *nphip,
int *sidehits
)
{
   float *phip,grid,lbound,hbound,phi,f;
   int nphi,i;
   struct range_list *rangep;
 
   setsdbound(srp,clp,&lbound,&hbound,&grid);
   if (!srp->vavoid)
   {
      /* ACRM added cast */
      nphi = *nphip = chmceil((f=(float)((hbound-lbound)/grid),&f));
      phip = *philistp = (float *) alloc((nphi+1) * sizeof(float));
      for (phi=lbound,i=1; phi<=hbound; phi += grid,i++)
      {
         if (i>nphi)
         {
            fprintf(out,"Error in setsidphis -- nphi exceeded.\n");
            die();
         }
         *phip++ = phi;
      }
   }
   else
   {
      explcntcts(srp,clp,&rangep,sidehits);
      mapgrdtorng(lbound,hbound,grid,rangep,philistp,nphip);
      free_range(rangep);
   }
}
 

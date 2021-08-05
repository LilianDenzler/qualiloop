#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Selects sidechain torsions which have at least one atom in the clumps
*   under srp.
*/
 
void selsidetor(
struct sideres *srp,
struct resphi **philist,
short ip[],
short jp[],
short kp[],
short lp[],
short icp[],
int nphi
)
{
   int iphi,count;
   struct resphi *rp;
   short int *iq,*jq,*kq,*lq,*icq;
 
   count = 0;
   for (iphi=0; iphi<nphi; iphi++)
   {
      walk_over_sideres(srp,
         if (atp->atomno == ip[iphi] ||
         atp->atomno == jp[iphi] ||
         atp->atomno == kp[iphi] ||
         atp->atomno == lp[iphi])
         {
            count++;
            goto nextphi1;
         }
      );
nextphi1:;
   }
   rp = *philist = (struct resphi *) alloc(sizeof(struct resphi));
   rp->nphi = count;
   iq = rp->ip = (short int *) alloc(count * sizeof(short int));
   jq = rp->jp = (short int *) alloc(count * sizeof(short int));
   kq = rp->kp = (short int *) alloc(count * sizeof(short int));
   lq = rp->lp = (short int *) alloc(count * sizeof(short int));
   icq = rp->icp = (short int *) alloc(count * sizeof(short int));
   for (iphi=0; iphi<nphi; iphi++)
   {
      walk_over_sideres(srp,
         if (atp->atomno == ip[iphi] ||
             atp->atomno == jp[iphi] ||
             atp->atomno == kp[iphi] ||
             atp->atomno == lp[iphi])
         {
            *iq++ = ip[iphi];
            *jq++ = jp[iphi];
            *kq++ = kp[iphi];
            *lq++ = lp[iphi];
            *icq++ = icp[iphi];
            goto nextphi2;
         }
      );
nextphi2:;
   }
}
 

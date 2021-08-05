#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Checks array boundaries for the residue topology information.
   02.08.92 Some recoding.    By:   ACRM
*/
 
void CheckRTF(void)
{
   int  i;
   char type[16];
   
   if(restop.nreses != 0)
   {
      for(i=0;i<ninfo;i++)
      {
         if(restop.nparam[restop.nreses-1][i] > restop.nparmx[i])
         {
            switch(i)
            {
            case 0:
               strcpy(type, "Atoms");
               break;
            case 1:
               strcpy(type, "Bonds");
               break;
            case 2:
               strcpy(type, "Angles");
               break;
            case 3:
               strcpy(type, "Torsions");
               break;
            case 4:
               strcpy(type, "Impropers");
               break;
            case 5:
               strcpy(type, "Exclusions");
               break;
            case 6:
               strcpy(type, "Donors");
               break;
            case 7:
               strcpy(type, "Acceptors");
               break;
            case 8:
               strcpy(type, "Builds");
               break;
            case 9:
               strcpy(type, "Groups");
               break;
            }

            fprintf(out,"Error==> While reading residue topology, array \
dimension for %s is too small.\n",type);
            fprintf(out,"Residue number is: %2d; Number of %s is: %5d; \
Maximum is: %5d\n", restop.nreses,type,restop.nparam[restop.nreses][i],
restop.nparmx[i]);

            die();
         }

         if(restop.nparam[restop.nreses-1][i] > restop.nparus[i])
            restop.nparus[i] = restop.nparam[restop.nreses-1][i];
      }
   }
}
 

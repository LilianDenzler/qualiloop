#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
#define MAXBUFF 200
#ifndef PI
#define PI (4.0 * atan((double)1.0))
#endif
#define RAD (PI/180.0)
 
extern int info_level;

/* This routine processes angle specifications.
   04.08.92 Some rewriting.   By:   ACRM.
*/
 
void procang(
FILE *fp,
short IndxTab[100][100]
)
{
   char  buff1[MAXBUFF+1],
         buff2[MAXBUFF+1],
         *buffp,
         Atom1[8],
         Atom2[8],
         Atom3[8];
   float ForceConst,
         OptAngle;
   int   i,j,k,l;
 
   if(info_level) fprintf(out,
"\n         BOND ANGLE PARAMETERS\n         =====================\n\n\
       #      ANGLE           CODE        KT        TO\n\n");
 
   while(fgets(buff1,MAXBUFF,fp))
   {
      terminate(buff1);
      struppr(buff1,buff2);
      buffp = killspcs(buff2);
      if(!strncmp(buffp,"END",3)) break;
 
      sscanf(buffp,"%s %s %s %f %f",Atom1,Atom2,Atom3,&ForceConst,&OptAngle);
      ljustpad(Atom1);
      ljustpad(Atom2);
      ljustpad(Atom3);
 
      i = j = k = -1;
      for(l=0;l<values.natyps;l++)
      {
         if(!strncmp(restop.acodes[l],Atom1,4)) i = l;
         if(!strncmp(restop.acodes[l],Atom2,4)) j = l;
         if(!strncmp(restop.acodes[l],Atom3,4)) k = l;
      }
      if(i == -1 || j == -1 || k == -1)
      {
         fprintf(out,"Warning==> Atoms in angle %4s %4s %4s %10.5f %10.5f \
do not exist in residue topology.\n",Atom1,Atom2,Atom3,ForceConst,OptAngle);
      }
      else
      {
         char temp1[5],
              temp2[5],
              temp3[5];
 
         engpar.angkey[values.napar] = values.natyps*(int)IndxTab[i][k]+j+1;
         engpar.angcon[values.napar] = ForceConst;
         engpar.eqang[values.napar] = OptAngle * RAD;
 
         strncpy(temp1,restop.acodes[i],4); temp1[4] = '\0';
         strncpy(temp2,restop.acodes[j],4); temp2[4] = '\0';
         strncpy(temp3,restop.acodes[k],4); temp3[4] = '\0';

         if(info_level) fprintf(out,"      %3d  %4s - %4s - %4s%9d%10.2f%10.2f\n",
                 values.napar+1,temp1,temp2,temp3,
                 engpar.angkey[values.napar],ForceConst,OptAngle);
         values.napar++;
      }
   }
}
 

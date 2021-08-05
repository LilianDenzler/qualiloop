#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
#define MAXBUFF 200

#ifndef PI
#define PI (4.0*atan(1.0))
#endif

#define RAD (PI/180.0)

#define ATOM_X "X   "
 
extern int info_level;

/* Process torsion angle information.
   02.08.92 Some rewriting.   By:   ACRM
*/
 
void proctors(
FILE *fp,
short IndxTab[100][100]
)
{
   char buff1[MAXBUFF+1],
        buff2[MAXBUFF+1],
        *buffp,
        atom_i[8],
        atom_j[8],
        atom_k[8],
        atom_l[8],
        atom1[8],
        atom4[8],
        temp1[8],
        temp2[8];
   int  i1,j1,k1,l1,j;
   float TorOptimum,Periodicity,ForceConst;
 
   if(info_level) fprintf(out,"\n                     DIHEDRAL ANGLE PARAMETERS\n\
                     =========================\
\n\n       #         TORSION               CODE        KP         N       \
DELTA\n");
 
   while(fgets(buff1,MAXBUFF,fp))
   {
      terminate(buff1);
      struppr(buff1,buff2);
      buffp = killspcs(buff2);
      if(!strncmp(buffp,"END",3)) break;
 
      sscanf(buffp,"%s %s %s %s %f %f %f",
             atom_i,atom_j,atom_k,atom_l,&ForceConst,&Periodicity,&TorOptimum);
      ljustpad(atom_i);
      ljustpad(atom_j);
      ljustpad(atom_k);
      ljustpad(atom_l);

      i1 = j1 = k1 = l1 = -1;
      for(j=0;j<values.natyps;j++)
      {
         if(!strncmp(restop.acodes[j],atom_i,4)) k1 = j;
         if(!strncmp(restop.acodes[j],atom_j,4)) i1 = j;
         if(!strncmp(restop.acodes[j],atom_k,4)) j1 = j;
         if(!strncmp(restop.acodes[j],atom_l,4)) l1 = j;
      }
 
      if(i1 >= 0 && j1 >= 0)
      {
         if(k1 == -1 && l1 == -1)
         {  /* Both unknown */
            strncpy(atom1,ATOM_X,4);      atom1[4] = '\0';
            strncpy(atom4,ATOM_X,4);      atom4[4] = '\0';

            engpar.torkey[values.nptpar] = (int)IndxTab[i1][j1];
            engpar.torcon[values.nptpar] = ForceConst;
            engpar.tormlt[values.nptpar] = Periodicity;
            engpar.torphs[values.nptpar] = TorOptimum*RAD;

            strncpy(temp1,restop.acodes[i1],4); temp1[4] = '\0';
            strncpy(temp2,restop.acodes[j1],4); temp2[4] = '\0';

            if(info_level) fprintf(out,"      %3d  %4s - %4s - %4s - %4s\
%9d%10.2f%10.2f%10.2f\n",values.nptpar+1,atom1,temp1,temp2,atom4,
engpar.torkey[values.nptpar],ForceConst,Periodicity,TorOptimum);
            values.nptpar++;
         }
         else if(k1 >= 0 && l1 >= 0)
         {  /* Both known */
            strncpy(atom1,restop.acodes[k1],4); atom1[4] = '\0';
            strncpy(atom4,restop.acodes[l1],4); atom4[4] = '\0';

            engpar.torkey[values.nptpar] = (int)IndxTab[i1][j1] +
               (int)IndxTab[k1][l1]*values.natyps*(values.natyps+1)/2;
            engpar.torcon[values.nptpar] = ForceConst;
            engpar.tormlt[values.nptpar] = Periodicity;
            engpar.torphs[values.nptpar] = TorOptimum*RAD;

            strncpy(temp1,restop.acodes[i1],4); temp1[4] = '\0';
            strncpy(temp2,restop.acodes[j1],4); temp2[4] = '\0';

            if(info_level) fprintf(out,"      %3d  %4s - %4s - %4s - %4s\
%9d%10.2f%10.2f%10.2f\n",values.nptpar+1,atom1,temp1,temp2,atom4,
engpar.torkey[values.nptpar],ForceConst,Periodicity,TorOptimum);
            values.nptpar++;
         }
         else
         {  /* One known, but not the other */
            fprintf(out,"Warning==> Atoms in PHI %4s - %4s - %4s - %4s\
%10.5f%10.5f%10.5f don't exist\n",
atom_i,atom_j,atom_k,atom_l,ForceConst,Periodicity,TorOptimum);
         }
      }
      else
      {
         fprintf(out,"Warning==> Atoms in PHI %4s - %4s - %4s - %4s\
%10.5f%10.5f%10.5f don't exist\n",
atom_i,atom_j,atom_k,atom_l,ForceConst,Periodicity,TorOptimum);
      }
   }
}
 

#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
#define MAXBUFF 200
 
extern int info_level;
 
void procimps(
FILE *fp,
short IndxTab[100][100]
)
{
   char buff1[MAXBUFF+1],
        buff2[MAXBUFF+1],
        *buffp,
        atom_i[5],
        atom_j[5],
        atom_k[5],
        atom_l[5],
        atck[8],
        atcl[8],
        ax[5],
        temp1[5],
        temp2[5];
   int  i1,j1,k1,l1,j;
   float tor0,fconst;
 
   strcpy(ax,"X   ");
 
   if(info_level) fprintf(out,"\n         I M P R O P E R   T O R S I O N   P A R A M E T E R S\
\n\n       #         T O R S I O N         CODE        KP         N\n");
 
   while(fgets(buff1,MAXBUFF,fp))
   {
      terminate(buff1);
      struppr(buff1,buff2);
      buffp = killspcs(buff2);
      if(!strncmp(buffp,"END",3)) break;
 
      sscanf(buffp,"%s %s %s %s %f %f",
             atom_i,atom_k,atom_l,atom_j,&fconst,&tor0);
      ljustpad(atom_i);
      ljustpad(atom_j);
      ljustpad(atom_k);
      ljustpad(atom_l);
      i1 = j1 = k1 = l1 = -1;
      for(j=0;j<values.natyps;j++)
      {
         if(!strncmp(restop.acodes[j],atom_i,4)) i1 = j;
         if(!strncmp(restop.acodes[j],atom_j,4)) j1 = j;
         if(!strncmp(restop.acodes[j],atom_k,4)) k1 = j;
         if(!strncmp(restop.acodes[j],atom_l,4)) l1 = j;
      }
 
      if(i1 >= 0 && j1 >= 0)
      {
         if(k1 == -1 && l1 == -1)
         {  /* Both unknown */
            strncpy(atcl,ax,4);
            atcl[4] = '\0';
            strncpy(atck,ax,4);
            atck[4] = '\0';
            engpar.impkey[values.nitpar] = (int)IndxTab[i1][j1];
            engpar.impcon[values.nitpar] = fconst;
            engpar.eqitan[values.nitpar] = tor0 * RAD;
            strncpy(temp1,restop.acodes[i1],4); temp1[4] = '\0';
            strncpy(temp2,restop.acodes[j1],4); temp2[4] = '\0';
            if(info_level) fprintf(out,"      %3d  %4s - %4s - %4s - %4s \
%9d %10.2f %10.2f\n",values.nitpar+1,temp1,atck,atcl,temp2,
engpar.impkey[values.nitpar],fconst,tor0);
            values.nitpar++;
         }
         else if(k1 >= 0 && l1 >= 0)
         {  /* Both known */
            strncpy(atcl,restop.acodes[l1],4);
            atcl[4] = '\0';
            strncpy(atck,restop.acodes[k1],4);
            atck[4] = '\0';
            engpar.impkey[values.nitpar] = (int)IndxTab[i1][j1] +
               (int)IndxTab[k1][l1]*values.natyps*(values.natyps+1)/2;
            engpar.impcon[values.nitpar] = fconst;
            engpar.eqitan[values.nitpar] = tor0 * RAD;
            strncpy(temp1,restop.acodes[i1],4); temp1[4] = '\0';
            strncpy(temp2,restop.acodes[j1],4); temp2[4] = '\0';
            if(info_level) fprintf(out,"      %3d  %4s - %4s - %4s - %4s \
%9d %10.2f %10.2f\n",values.nitpar+1,temp1,atck,atcl,temp2,
engpar.impkey[values.nitpar],fconst,tor0);
            values.nitpar++;
         }
         else
         {  /* One known, but not the other */
            printf("Warning: Atoms in IMPHI %4s - %4s - %4s - %4s \
%10.5f %10.5f don't exist\n",
atom_i,atom_j,atom_k,atom_l,fconst,tor0);
         }
      }
      else
      {
         printf("Warning: Atoms in IMPHI %4s - %4s - %4s - %4s \
%10.5f %10.5f don't exist\n",
atom_i,atom_j,atom_k,atom_l,fconst,tor0);
      }
   }
}
 

#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
#define MAXBUFF 200
 
int ReadParams(
FILE *fp
)
{
   char buff1[MAXBUFF+1],
        buff2[MAXBUFF+1],
        *buffp;
   short IndxTab[100][100];
   int i,j,k,
       ncn,
       idx[500];
   float temp1[500],
         temp2[500],
         temp3[500];
 
   /* Set up IndxTab[][], the temporary code numbers array */
   k = 0;
   for(i=1;i<=values.natyps;i++)
      for(j=1;j<=i;j++)
      {
         k++;
         IndxTab[i-1][j-1] = (short)k;
         IndxTab[j-1][i-1] = (short)k;
      }
 
   /* Set various starting values */
   values.nbpar = 0;
   values.napar = 0;
   values.nptpar = 0;
   values.nitpar = 0;
   values.nhbpar = 0;
 
   ncn = 0;
 
   rdprtitle(fp);
   while(fgets(buff1,MAXBUFF,fp))
   {
      terminate(buff1);
      struppr(buff1,buff2);
      buffp = killspcs(buff2);
 
      if(!strncmp(buffp,"BOND",4))
         ProcessBonds(fp,IndxTab);
      if(!strncmp(buffp,"ANGL",4))
         procang(fp,IndxTab);
      if(!strncmp(buffp,"TORS",4))
         proctors(fp,IndxTab);
      if(!strncmp(buffp,"IMPR",4))
         procimps(fp,IndxTab);
      if(!strncmp(buffp,"NONB",4))
         procnbnds(fp,IndxTab,&ncn);
      if(!strncmp(buffp,"HBON",4))
         prochbnd(fp,IndxTab);
   }
   /* Now sort the code arrays by key number */
   sorti_perm(engpar.bndkey,values.nbpar,idx);
   for(i=0;i<values.nbpar;i++)
   {
      temp1[i] = engpar.bndcon[idx[i]];
      temp2[i] = engpar.eqbdis[idx[i]];
   }
   for(i=0;i<values.nbpar;i++)
   {
      engpar.bndcon[i] = temp1[i];
      engpar.eqbdis[i] = temp2[i];
   }
   sorti_perm(engpar.angkey,values.napar,idx);
   for(i=0;i<values.napar;i++)
   {
      temp1[i] = engpar.angcon[idx[i]];
      temp2[i] = engpar.eqang[idx[i]];
   }
   for(i=0;i<values.napar;i++)
   {
      engpar.angcon[i] = temp1[i];
      engpar.eqang[i] = temp2[i];
   }
 
   sorti_perm(engpar.torkey,values.nptpar,idx);
   for(i=0;i<values.nptpar;i++)
   {
      temp1[i] = engpar.torcon[idx[i]];
      temp2[i] = engpar.tormlt[idx[i]];
      temp3[i] = engpar.torphs[idx[i]];
   }
   for(i=0;i<values.nptpar;i++)
   {
      engpar.torcon[i] = temp1[i];
      engpar.tormlt[i] = temp2[i];
      engpar.torphs[i] = temp3[i];
   }
 
   sorti_perm(engpar.impkey,values.nitpar,idx);
   for(i=0;i<values.nitpar;i++)
   {
      temp1[i] = engpar.impcon[idx[i]];
      temp2[i] = engpar.eqitan[idx[i]];
   }
   for(i=0;i<values.nitpar;i++)
   {
      engpar.impcon[i] = temp1[i];
      engpar.eqitan[i] = temp2[i];
   }
 
   sorti_perm(engpar.hbkey,values.nhbpar,idx);
   for(i=0;i<values.nhbpar;i++)
   {
      temp1[i] = engpar.hbr12[idx[i]];
      temp2[i] = engpar.hbr10[idx[i]];
   }
   for(i=0;i<values.nhbpar;i++)
   {
      engpar.hbr12[i] = temp1[i];
      engpar.hbr10[i] = temp2[i];
   }
 
   /* We won't need the file again, so close it */
   fclose(fp);
   fp = NULL;
 
   return(ncn);
}
 

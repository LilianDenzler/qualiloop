/* This routine sets up the codes table for rapid access to bond,
   angle, torsion, improper, and hbond parameters.
   
   02.08.92 Recoded     By:   ACRM
*/

#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif
 
void codes(int DoPstruct,
           int DoHBonds)
{
   int   PhiFlag       = TRUE,
         ImproperFlag  = TRUE,
         DieFlag       = FALSE,
         IndxFactor,
         OldICP        = 0,
         i,    j,    k = 0,
         i1,   j1,   k1,   l1,
         IndxPos,
         IndxPos2,
         NumAtTypeCodes = values.natyps,
         NumBonds       = values.nbonds,
         NumBondCodes   = values.nbpar,
         NumAngles      = values.nangs,
         NumAngleCodes  = values.napar,
         NumTor         = values.nptors,
         NumTorCodes    = values.nptpar,
         NumImp         = values.nitors,
         NumImpCodes    = values.nitpar,
         NumHBond       = values.nhbs,
         NumHBondCodes  = values.nhbpar;

   char  temp1[8],
         temp2[8],
         temp3[8];
 
   short IndxTab[100][100];
 
   IndxFactor = NumAtTypeCodes*(NumAtTypeCodes+1)/2;
 
   for(i=0; i<NumAtTypeCodes; i++)
      for(j=0; j<=i; j++)
         IndxTab[i][j] = IndxTab[j][i] = (int)(++k);
 
/*******************************
 *** Deal with PSTRUCT stuff ***
 *******************************/
   if(DoPstruct)
   {
/*** Bonds ***/
      if(NumBonds)
      {
         for(i=0; i<NumBonds; i++)
         {
            if(pstruct.atbnd1[i] > 0)
            {
               i1   = pstruct.atcode[pstruct.atbnd1[i]-1] - 1;
               j1   = pstruct.atcode[pstruct.atbnd2[i]-1] - 1;
               IndxPos = (int)IndxTab[j1][i1];

               if(!(getpar.ibndp[i] = bin_search(IndxPos,engpar.bndkey,NumBondCodes)))
               {
                  strncpy(temp1, restop.acodes[i1], 4); temp1[4] = '\0';
                  strncpy(temp2, restop.acodes[j1], 4); temp2[4] = '\0';
                  fprintf(out, "Error==> Parameters not found for bond %d \
between atoms %7s %7s\n         Program dies!\n", i+1, temp1, temp2);
                  DieFlag = TRUE;
               }
            }
         }
      }
 
/*** Angles ***/
      if(NumAngles)
      {
         for(i=0; i<NumAngles; i++)
         {
            if(pstruct.atang1[i] > 0)
            {
               i1   = pstruct.atcode[pstruct.atang1[i]-1] - 1;
               j1   = pstruct.atcode[pstruct.atang2[i]-1] - 1;
               k1   = pstruct.atcode[pstruct.atang3[i]-1] - 1;
               IndxPos = NumAtTypeCodes * (int)IndxTab[k1][i1] + j1 + 1;

               if(!(getpar.iangp[i] = bin_search(IndxPos,engpar.angkey,
NumAngleCodes)))
               {
                  strncpy(temp1, restop.acodes[i1], 4); temp1[4] = '\0';
                  strncpy(temp2, restop.acodes[j1], 4); temp2[4] = '\0';
                  strncpy(temp3, restop.acodes[k1], 4); temp3[4] = '\0';
                  fprintf(out, "Error==> Parameters not found for angle %d \
between atoms %7s %7s %7s\n         Program dies!\n", i+1, temp1, temp2, temp3);
                  DieFlag = TRUE;
               }
            }
         }
      }
 
/*** Proper torsions ***/
      if(NumTor)
      {
         for(i=0; i<NumTor; i++)
         {
            if(pstruct.attor1[i] > 0)
            {
               i1    = pstruct.atcode[pstruct.attor1[i]-1] - 1;
               j1    = pstruct.atcode[pstruct.attor2[i]-1] - 1;
               k1    = pstruct.atcode[pstruct.attor3[i]-1] - 1;
               l1    = pstruct.atcode[pstruct.attor4[i]-1] - 1;
               IndxPos  = (int)IndxTab[k1][j1];
               IndxPos2 = IndxPos+(int)IndxTab[l1][i1]*IndxFactor;

               if((getpar.itorp[i] = bin_search(IndxPos2, engpar.torkey,
 NumTorCodes)))
               {
                  getpar.itorp[i] = 0;
                  for(j=0; j<NumTorCodes; j++)
                  {
                     if(IndxPos2 == engpar.torkey[j] && j+1 != OldICP)
                        getpar.itorp[i] = j+1;
                  }

                  OldICP=0;
                  if((i+1 != NumTor) &&
                     (pstruct.attor1[i]==pstruct.attor1[i+1]) &&
                     (pstruct.attor2[i]==pstruct.attor2[i+1]) &&
                     (pstruct.attor3[i]==pstruct.attor3[i+1]) &&
                     (pstruct.attor4[i]==pstruct.attor4[i+1]))
                  {
                     if(PhiFlag)
                     {
                        fprintf(out, "Note==> Using multiple potential terms \
for some dihedrals\n");
                        PhiFlag = FALSE;
                     }
                     OldICP = getpar.itorp[i];
                  }
               }
               else
               {
                  OldICP = 0;
                  getpar.itorp[i] = bin_search(IndxPos, engpar.torkey,
 NumTorCodes);
               }

               if(getpar.itorp[i] == 0)
               {
                  strncpy(temp1, restop.acodes[j1], 4); temp1[4] = '\0';
                  strncpy(temp2, restop.acodes[k1], 4); temp2[4] = '\0';
                  fprintf(out, "Error==> Parameters not found for torsion %d \
between atoms %7s %7s\n         Program dies!\n", i+1, temp1, temp2);
                  DieFlag = TRUE;
               }
            }
         }
         OldICP = 0;
      }
 
/*** Impropers ***/
      if(NumImp)
      {
         for(i=0; i<NumImp; i++)
         {
            if(pstruct.atimp1[i] > 0)
            {
               i1    = pstruct.atcode[pstruct.atimp1[i]-1]-1;
               j1    = pstruct.atcode[pstruct.atimp2[i]-1]-1;
               k1    = pstruct.atcode[pstruct.atimp3[i]-1]-1;
               l1    = pstruct.atcode[pstruct.atimp4[i]-1]-1;
               IndxPos  = (int)IndxTab[l1][i1];
               IndxPos2 = IndxPos+(int)(IndxTab[k1][j1])*IndxFactor;

               if((getpar.iimpp[i] = bin_search(IndxPos2, engpar.impkey,
 NumImpCodes)))
               {
                  getpar.iimpp[i]=0;
                  for(j=0; j<NumImpCodes; j++)
                  {
                     if(IndxPos2 == engpar.impkey[j] &&
                        j+1 != OldICP) getpar.iimpp[i] = j+1;
                  }
                  
                  OldICP = 0;
                  if((i+1 != NumImp) &&
                     (pstruct.atimp1[i]==pstruct.atimp1[i+1]) &&
                     (pstruct.atimp2[i]==pstruct.atimp2[i+1]) &&
                     (pstruct.atimp3[i]==pstruct.atimp3[i+1]) &&
                     (pstruct.atimp4[i]==pstruct.atimp4[i+1]))
                  {
                     if(ImproperFlag)
                     {
                        fprintf(out, "Note==> Using multiple potential terms \
for some impropers\n");
                        ImproperFlag = FALSE;
                     }
                     OldICP = getpar.iimpp[i];
                  }
               }
               else
               {
                  OldICP = 0;
                  getpar.iimpp[i] = bin_search(IndxPos, engpar.impkey,
 NumImpCodes);
               }

               if(getpar.iimpp[i] == 0)
               {
                  strncpy(temp1, restop.acodes[i1], 4); temp1[4] = '\0';
                  strncpy(temp2, restop.acodes[l1], 4); temp2[4] = '\0';
                  fprintf(out, "Error==> Parameters not found for improper %d \
between atoms %7s %7s\n         Program dies!\n", i+1, temp1, temp2);
                  DieFlag = TRUE;
               }
            }
         }
      }
   }
   
/*************************
 *** Deal with H bonds ***
 *************************/
   if(DoHBonds && NumHBond)
   {
      for(i=0; i<NumHBond; i++)
      {
         if(hbonds.hbdon[i] > 0)
         {
            i1   = pstruct.atcode[hbonds.hbdon[i]-1] - 1;
            j1   = pstruct.atcode[hbonds.hbacc[i]-1] - 1;
            IndxPos = (int)IndxTab[j1][i1];

            if(!(getpar.ihbp[i] = bin_search(IndxPos,engpar.hbkey,
NumHBondCodes)))
            {
               strncpy(temp1, restop.acodes[i1], 4); temp1[4] = '\0';
               strncpy(temp2, restop.acodes[j1], 4); temp2[4] = '\0';
               fprintf(out, "Error==> Parameters not found for hydrogen bond %d \
between atoms %7s %7s\n         Program dies!\n", i+1, temp1, temp2);
               DieFlag = TRUE;
            }
         }
      }
 
   }

   if(DieFlag) die();
}
 

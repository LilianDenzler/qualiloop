#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

/* Generates the internal coordinate lists
   02.08.92 Some recoding.    By:   ACRM.
*/

void genic(int Start,
           int Stop)
{
   char temp[8];
   int i, j, k, l, iblold;
   int NumAtoms, 
       NumBonds, 
       NumAngles, 
       NumTorsions, 
       NumImpropers, 
       NumNBonds, 
       NumDonors, 
       NumAcceptors, 
       nl, itemp;
 
   /* Search the residue topology information for this residue 
      from the protein structure information.
   */
   for(i=Start; i<=Stop; i++)
   {
      int found = FALSE;
 
      strncpy(temp,pstruct.resnme[i-1],4);
      for(j=1; j<=restop.nreses; j++)
      {
         if(!strncmp(temp,restop.namres[j-1],4))
         {
            pstruct.resndx[i-1] = j;
            found = TRUE;
            break;
         }
      }
      
      if(!found)
      {
         temp[4] = '\0';
         fprintf(out,"Error==> Unrecognised Residue %4s while generating \
internal coordinates.\n",temp);
         die();
      }
   }
 
   /* Generate the internal coordinates */
   iblold = 0;
   if(values.natoms > 0) iblold = pstruct.lstexc[values.natoms-1];
   for(j=Start; j<=Stop; j++)
   {
      pstruct.lstatm[j-1] = values.natoms;
      k              = pstruct.resndx[j-1];
      NumAtoms       = restop.nparam[k-1][0];
      NumBonds       = restop.nparam[k-1][1];
      NumAngles      = restop.nparam[k-1][2];
      NumTorsions    = restop.nparam[k-1][3];
      NumImpropers   = restop.nparam[k-1][4];
      NumNBonds      = restop.nparam[k-1][5];
      NumDonors      = restop.nparam[k-1][6];
      NumAcceptors   = restop.nparam[k-1][7];
 
      if(NumBonds)
      {
         for(l=0; l<NumBonds; l++)
         {
            values.nbonds++;
            pstruct.atbnd1[values.nbonds-1] = values.natoms + restop.bndat1[k-1][l];
            pstruct.atbnd2[values.nbonds-1] = values.natoms + restop.bndat2[k-1][l];
         }
      }
      
      if(NumAngles)
      {
         for(l=0; l<NumAngles; l++)
         {
            values.nangs++;
            pstruct.atang1[values.nangs-1] = values.natoms+restop.angat1[k-1][l];
            pstruct.atang2[values.nangs-1] = values.natoms+restop.angat2[k-1][l];
            pstruct.atang3[values.nangs-1] = values.natoms+restop.angat3[k-1][l];
         }
      }

      if(NumTorsions)
      {
         for(l=0; l<NumTorsions; l++)
         {
           values.nptors++;
           pstruct.attor1[values.nptors-1] = values.natoms+restop.torat1[k-1][l];
           pstruct.attor2[values.nptors-1] = values.natoms+restop.torat2[k-1][l];
           pstruct.attor3[values.nptors-1] = values.natoms+restop.torat3[k-1][l];
           pstruct.attor4[values.nptors-1] = values.natoms+restop.torat4[k-1][l];
         }
      }

      if(NumImpropers)
      {
         for(l=0; l<NumImpropers; l++)
         {
            values.nitors++;
           pstruct.atimp1[values.nitors-1] = values.natoms+restop.impat1[k-1][l];
           pstruct.atimp2[values.nitors-1] = values.natoms+restop.impat2[k-1][l];
           pstruct.atimp3[values.nitors-1] = values.natoms+restop.impat3[k-1][l];
           pstruct.atimp4[values.nitors-1] = values.natoms+restop.impat4[k-1][l];
         }
      }

      if(NumNBonds)
      {
         for(l=0; l<NumNBonds; l++)
         {
           values.nnbs++;
           if(restop.exclnb[k-1][l] == -99)
             pstruct.nbexcl[values.nnbs-1] = 0;
           else
             pstruct.nbexcl[values.nnbs-1] = values.natoms+restop.exclnb[k-1][l];
         }
      }

      if(NumAtoms)
      {
         for(l=0; l<NumAtoms; l++)
         {
            nl                 = values.natoms+l;
            pstruct.lstexc[nl] = restop.nexcld[k-1][l] + iblold;
            iblold             = pstruct.lstexc[nl];
            itemp              = restop.acindx[k-1][l];
            pstruct.atcode[nl] = itemp;
            pstruct.atmass[nl] = restop.atmmas[itemp-1];
            pstruct.atchrg[nl] = restop.atmchg[k-1][l];
            strncpy(pstruct.atmnme[nl],restop.namatm[k-1][l],4);
         }
      }

      if(NumAcceptors)
      {
         for(l=0; l<NumAcceptors; l++)
         {
            values.naccat++;

            if(restop.aan1hb[k-1][l] == -99)
               pstruct.hbaan1[values.naccat-1] = 0;
            else
               pstruct.hbaan1[values.naccat-1] =
                            values.natoms + restop.aan1hb[k-1][l];

            if(restop.aan2hb[k-1][l] == -99)
               pstruct.hbaan2[values.naccat-1] = 0;
            else
               pstruct.hbaan2[values.naccat-1] =
                            values.natoms + restop.aan2hb[k-1][l];

            pstruct.hbacpt[values.naccat-1] = 
                         values.natoms + restop.acpthb[k-1][l];
         }
      }

      if(NumDonors)
      {
         for(l=0; l<NumDonors; l++)
         {
            values.ndonat++;
            
            if(restop.dhydhb[k-1][l] == -99)
               pstruct.hbdhyd[values.ndonat-1] = 0;
            else
               pstruct.hbdhyd[values.ndonat-1] =
                            values.natoms + restop.dhydhb[k-1][l];

            if(restop.dan1hb[k-1][l] == -99)
               pstruct.hbdan1[values.ndonat-1] = 0;
            else
               pstruct.hbdan1[values.ndonat-1] =
                            values.natoms + restop.dan1hb[k-1][l];

            if(restop.dan2hb[k-1][l] == -99)
               pstruct.hbdan2[values.ndonat-1] = 0;
            else
               pstruct.hbdan2[values.ndonat-1] = 
                            values.natoms + restop.dan2hb[k-1][l];

            pstruct.hbdonr[values.ndonat-1] =
                         values.natoms + restop.donrhb[k-1][l];
         }
         /* Original code puts this here, I don't believe it! */
         /* values.natom += NumAtoms; */
      }
      values.natoms += NumAtoms;
   }
   pstruct.lstatm[Stop] = values.natoms;
}
 

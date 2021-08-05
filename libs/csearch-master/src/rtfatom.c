#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Reads information for an atom. The exclusion names must be in double
   inverted commas if more than one.
 
   Syntax: ATOM iupac atom-type-code charge "repeat(exclusion-names)"

   04.08.92 Some rewriting.   By:   ACRM.
*/
 
void RTFAtom(
int  *NBondExclNum,
int  SawTable,
int  *NumNBondExcl,
char NBondExclNam[mxnber][4],
int  *NumAtTypes
)
{
 
   int  NumAtoms,
        i, 
        ok=TRUE, 
        TypeCode, 
        OldNumNBExcl;
   char *buffer, 
        name[8];
 
   abmpad(ResTop_strparam[0],4);
   abmpad(ResTop_strparam[1],4);
 
   if(restop.nreses == 0)
   {
      fprintf(out,"Error==> A residue must be specified before atoms in the \
residue topology.\n");
      ResTop_errcnt++;
   }
   else
   {
      if(SawTable)
      {
         fprintf(out,"Error==> An atom may not be defined in the residue topology\n");
         fprintf(out,"         once information in the DECLARE table is used in a\n");
         fprintf(out,"         residue as it will change the number of atoms in a\n");
         fprintf(out,"         residue. The atom will be ignored.\n");
         ResTop_errcnt++;
      }
      else
      {
         NumAtoms = restop.nparam[restop.nreses-1][0] + 1;
         if(!strncmp(ResTop_strparam[0],BLANK,4))
         {
            fprintf(out,"Error==> Atoms must have names in residue topology.\n");
            ResTop_errcnt++;
         }
         else
         {
            /* Check for dupes */
            for(i=0;i<NumAtoms-1;i++)
            {
               if(!strncmp(restop.namatm[restop.nreses-1][i],ResTop_strparam[0],4))
               {
                  fprintf(out,"Duplicate IUPAC atom name\n");
                  ResTop_errcnt++;
                  ok = FALSE;
                  break;
               }
            }

            if(ok)
            {
               /* Now we can actually set the atom name             */
               strncpy(restop.namatm[restop.nreses-1][NumAtoms-1],ResTop_strparam[0],4);

               if(!strncmp(ResTop_strparam[1],BLANK,4))
                  TypeCode = 0;
               else
               {
                  TypeCode = 0;
                  for(i=0; i<*NumAtTypes; i++)
                  {
                     if(!strncmp(restop.acodes[i],ResTop_strparam[1],4))
                     {
                        TypeCode = i+1;
                        break;
                     }
                  }
               }

               if(!TypeCode)
               {
                  fprintf(out,"Error==> Unknown atom type code in residue \
topology. Atom ignored.\n");
                  ResTop_errcnt++;
               }
               else
               {
                  int NumExcl = 0;
 
                  restop.acindx[restop.nreses-1][NumAtoms-1] = TypeCode;
                  sscanf(ResTop_strparam[2],"%f",&restop.atmchg[restop.nreses-1][NumAtoms-1]);
                  OldNumNBExcl = *NumNBondExcl;
                  buffer       = ResTop_strparam[3];

                  for(;;)
                  {
                     buffer += GetString(buffer,name);
                     buffer  = killspcs(buffer);
 
/*                   if(!name[0]) break; */
                     if(NumExcl && !name[0]) break;
 
                     abmpad(name,4);
                     if(strncmp(name,BLANK,4) || *NumNBondExcl==OldNumNBExcl)
                     {
                        /* Add a blank exclusion to keep GENIC happy.
                           Should be deleted when the routine which generates
                           internal coordinates can handle the condition of
                           no exclusions correctly.
                        */
                        if(++(*NumNBondExcl) >= mxnber)
                        {
                           fprintf(out,"Error==> Residue topology has too many \
non-bonded exclusions in one residue.\n");
                           fprintf(out,"         MXRTX (%d) exceeded.\n",
                                   mxnber);
                           die();
                        }
                        strncpy(NBondExclNam[*NumNBondExcl-1],name,4);
                        NumExcl++;
                     }
                  }
                  NBondExclNum[NumAtoms-1] = *NumNBondExcl;
                  restop.nparam[restop.nreses-1][0]  = NumAtoms;
               }
            }
         }
      }
   }
}
 

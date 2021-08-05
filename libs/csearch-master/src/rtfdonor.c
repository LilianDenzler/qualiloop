#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Specifies an HBond donor. This code only supports the explicit
   hydrogen model, so the first atom must be a hydrogen and
   antecedants must be specified.
   Syntax: DONOR hydrogen heavy-atom antecedant1 antecedant2

   04.08.92 Some rewriting.   By:   ACRM.
*/
 
void RTFDonor(
int  *p_SawTable,
int  TableSize,
char TableKey[MXTABL][4],
int  *TableReplace
)
{
   int NumDonor,
       Atom1, 
       Atom2, 
       Atom3, 
       Atom4, 
       ok, 
       NumAtoms, 
       HaveH = FALSE;

   if(restop.nreses == 0)
   {
      fprintf(out,"Error==> A residue must be specified before \
donors in the residue topology.\n");
      ResTop_errcnt++;
   }
   else
   {
      abmpad(ResTop_strparam[0],4);
      abmpad(ResTop_strparam[1],4);
      abmpad(ResTop_strparam[2],4);
      abmpad(ResTop_strparam[3],4);
 
      NumDonor = restop.nparam[restop.nreses-1][6];
      Atom1 = LookupName(ResTop_strparam[0],&ok,TableSize,TableKey,
                         TableReplace,p_SawTable);
      if(ok)
      {
         NumAtoms = restop.nparam[restop.nreses-1][0];
         if(Atom1<1 || Atom1>NumAtoms)
         {
            fprintf(out,"Error==> The HBond donor atoms in the residue \
topology must be in the current residue.\n         The donor will be ignored.\n");
            ResTop_errcnt++;
         }
         else
         {
            HaveH = (restop.atmmas[restop.acindx[restop.nreses-1][Atom1-1]-1] < 3.5) ?
                    TRUE : FALSE;
            if(!HaveH)
            {
               fprintf(out,"Error==> Only explicit H types in residue \
topology supported.\n");
               die();
            }
 
            Atom2 = LookupName(ResTop_strparam[1],&ok,TableSize,TableKey,
                               TableReplace,p_SawTable);
            if(ok)
            {
               Atom3 = LookupName(ResTop_strparam[2],&ok,TableSize,TableKey,
                                  TableReplace,p_SawTable);
               if(ok)
               {
                  Atom4 = LookupName(ResTop_strparam[3],&ok,TableSize,TableKey,
                                     TableReplace,p_SawTable);
                  if(ok)
                  {
                     NumDonor++;
                     restop.dhydhb[restop.nreses-1][NumDonor-1]  = Atom1;
                     restop.donrhb[restop.nreses-1][NumDonor-1]  = Atom2;
                     restop.dan1hb[restop.nreses-1][NumDonor-1] = Atom3;
                     restop.dan2hb[restop.nreses-1][NumDonor-1] = Atom4;
                     restop.nparam[restop.nreses-1][6]          = NumDonor;
                  }
               }
            }
         }
      }
   }
}
 

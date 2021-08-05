#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Specifies a bond between 2 atoms
   Syntax: BOND iupac1 iupac2

   04.08.92 Some rewriting.   By:   ACRM.
*/
 
void RTFBond(
int  *p_SawTable,
int  TableSize,
char TableKey[MXTABL][4],
int  TableReplace[MXTABL]
)
{
   int NumBonds,
       Atom1, Atom2,
       ok;

   if(restop.nreses == 0)
   {
      fprintf(out,"Error==> A residue must be specified in the residue topology \
before bonds.\n");
      ResTop_errcnt++;
   }
   else
   {
      NumBonds = restop.nparam[restop.nreses-1][1];
      abmpad(ResTop_strparam[0],4);
      abmpad(ResTop_strparam[1],4);
 
      Atom1 = LookupName(ResTop_strparam[0],&ok,TableSize,TableKey,
                         TableReplace,p_SawTable);
      if(ok)
      {
         Atom2 = LookupName(ResTop_strparam[1],&ok,TableSize,TableKey,
                            TableReplace,p_SawTable);
         if(ok)
         {
            if(Atom1 == Atom2)
            {
               fprintf(out,"Error==> Two identical atoms in residue topology \
bond will be ignored.\n");
               ResTop_errcnt++;
            }
            else
            {
               NumBonds++;
               restop.bndat1[restop.nreses-1][NumBonds-1] = Atom1;
               restop.bndat2[restop.nreses-1][NumBonds-1] = Atom2;
            }
         }
      }
      restop.nparam[restop.nreses-1][1] = NumBonds;
   }
}
 

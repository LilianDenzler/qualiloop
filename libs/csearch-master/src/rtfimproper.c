#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Specifies an improper torsion between 4 atoms
   Syntax: IMPROPER iupac1 iupac2 iupac3 iupac4

   04.08.92 Some rewriting.   By:   ACRM.
*/
 
void RTFImproper(
int  *p_SawTable,
int  TableSize,
char TableKey[MXTABL][4],
int  *TableReplace
)
{
   int NumImproper,
       Atom1, 
       Atom2, 
       Atom3, 
       Atom4, 
       ok;

   if(restop.nreses == 0)
   {
      fprintf(out,"Error==> A residue must be specified before \
impropers in residue topology.\n");
      ResTop_errcnt++;
   }
   else
   {
      abmpad(ResTop_strparam[0],4);
      abmpad(ResTop_strparam[1],4);
      abmpad(ResTop_strparam[2],4);
      abmpad(ResTop_strparam[3],4);
 
      NumImproper = restop.nparam[restop.nreses-1][4];
 
      Atom1 = LookupName(ResTop_strparam[0],&ok,TableSize,TableKey,
                         TableReplace,p_SawTable);
      if(ok)
      {
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
                  if(Atom1==Atom2 || Atom1==Atom3 || Atom1==Atom4 ||
                     Atom2==Atom3 || Atom2==Atom4 || Atom3==Atom4)
                  {
                     fprintf(out,"Error==> Atoms in an improper torsion in the \
residue topology are not all different.\n         The improper will be ignored.\n");
                     ResTop_errcnt++;
                  }
                  else
                  {
                     NumImproper++;
                     restop.impat1[restop.nreses-1][NumImproper-1] = Atom1;
                     restop.impat2[restop.nreses-1][NumImproper-1] = Atom2;
                     restop.impat3[restop.nreses-1][NumImproper-1] = Atom3;
                     restop.impat4[restop.nreses-1][NumImproper-1] = Atom4;
                  }
               }
            }
         }
      }
      restop.nparam[restop.nreses-1][4] = NumImproper;
   }
}
 

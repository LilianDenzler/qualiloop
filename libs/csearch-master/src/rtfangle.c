#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Specifies an angle between 3 atoms
   Syntax: ANGLE iupac1 iupac2 iupac3

   04.08.92 Some rewriting.   By:   ACRM.
*/
 
void RTFAngle(
int  *p_SawTable,
int  TableSize,
char TableKey[MXTABL][4],
int  *TableReplace
)
{
   int NumAngle, 
       Atom1, 
       Atom2, 
       Atom3, 
       ok;

   if(restop.nreses == 0)
   {
      fprintf(out,"Error==> A residue must be specified before bonds \
in residue topology.\n");
      ResTop_errcnt++;
   }
   else
   {
      abmpad(ResTop_strparam[0],4);
      abmpad(ResTop_strparam[1],4);
      abmpad(ResTop_strparam[2],4);
 
      NumAngle = restop.nparam[restop.nreses-1][2];
 
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
               if(Atom1==Atom2 || Atom1==Atom3 || Atom2==Atom3)
               {
                  fprintf(out,"Error==> Atoms in a residue topology angle are \
not all different.\n         The angle will be ignored.\n");
                  ResTop_errcnt++;
               }
               else
               {
                  NumAngle++;
                  restop.angat1[restop.nreses-1][NumAngle-1] = Atom1;
                  restop.angat2[restop.nreses-1][NumAngle-1] = Atom2;
                  restop.angat3[restop.nreses-1][NumAngle-1] = Atom3;
               }
            }
         }
      }
      restop.nparam[restop.nreses-1][2] = NumAngle;
   }
}
 

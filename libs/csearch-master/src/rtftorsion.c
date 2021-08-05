#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Specifies a torsion between 4 atoms
 
   Syntax: TORSION iupac1 iupac2 iupac3 iupac4

   02.08.92 Some Recoding.    By:   ACRM.
*/
 
void RTFTorsion(int  *p_SawTable,
                int  TableSize,
                char TableKey[MXTABL][4],
                int  *TableReplace)
{
   int NumTorsion, 
       TorPt1, TorPt2, 
       TorPt3, TorPt4,
       ok;

   if(restop.nreses == 0)
   {
      fprintf(out,"Error==> A residue must be specified before \
torsions in the residue topology.\n");
      ResTop_errcnt++;
   }
   else
   {
      abmpad(ResTop_strparam[0],4);
      abmpad(ResTop_strparam[1],4);
      abmpad(ResTop_strparam[2],4);
      abmpad(ResTop_strparam[3],4);
 
      NumTorsion = restop.nparam[restop.nreses-1][3];
 
      TorPt1 = LookupName(ResTop_strparam[0],&ok,TableSize,TableKey,
                          TableReplace,p_SawTable);
      if(ok)
      {
         TorPt2 = LookupName(ResTop_strparam[1],&ok,TableSize,TableKey,
                             TableReplace,p_SawTable);
         if(ok)
         {
            TorPt3 = LookupName(ResTop_strparam[2],&ok,TableSize,TableKey,
                                TableReplace,p_SawTable);
            if(ok)
            {
               TorPt4 = LookupName(ResTop_strparam[3],&ok,TableSize,TableKey,
                                   TableReplace,p_SawTable);
               if(ok)
               {
                  if(TorPt1==TorPt2 || TorPt1==TorPt3 || TorPt1==TorPt4 ||
                     TorPt2==TorPt3 || TorPt2==TorPt4 || TorPt3==TorPt4)
                  {
                     fprintf(out,"Error==> Atoms in a residue topology \
torsion are not all different.\n         The torsion will be ignored.\n");
                     ResTop_errcnt++;
                  }
                  else
                  {
                     NumTorsion++;
                     restop.torat1[restop.nreses-1][NumTorsion-1] = TorPt1;
                     restop.torat2[restop.nreses-1][NumTorsion-1] = TorPt2;
                     restop.torat3[restop.nreses-1][NumTorsion-1] = TorPt3;
                     restop.torat4[restop.nreses-1][NumTorsion-1] = TorPt4;
                  }
               }
            }
         }
      }
      restop.nparam[restop.nreses-1][3] = NumTorsion;
   }
}
 

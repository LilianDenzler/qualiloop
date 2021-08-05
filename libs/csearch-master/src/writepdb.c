#include "ProtoTypes.h"
#include "CongenProto.h"
 
void WritePDB(
FILE *fp,
PDB *pdb
)
{
   PDB *p;
 
   for(p = pdb ; p ; p = p->next)
   {
      wrtpdbrec(fp,p);
   }
}
 
 

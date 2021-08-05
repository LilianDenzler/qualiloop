#include "ProtoTypes.h"
#include "CongenProto.h"
 
void wrtpdbrec(
FILE *fp,
PDB *pdb
)
{
   fprintf(fp,"%6s%5d  %4s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
      pdb->junk,pdb->atnum,pdb->atnam,pdb->resnam,pdb->chain,pdb->resnum,
      pdb->insert,pdb->x,pdb->y,pdb->z,pdb->occ,pdb->bval);
}

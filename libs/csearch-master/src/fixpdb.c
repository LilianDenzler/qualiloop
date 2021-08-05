#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* This routine adds hydrogens to the PDB linked list and fixes the
   coords of CTER OT2 */
 
void fixpdb(
FILE *hadd_fp,
PDB *pdb
)
{
   PDB *p;
   int atnum = 0;
 
   hadd(hadd_fp,pdb);
 
   FixNterH(pdb);
   FixCterO(pdb);
 
   /* Renumber the pdb linked list */
   for(p=pdb; p; NEXT(p)) p->atnum = ++atnum;
 }
 

#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"

/* Read header info from the CG file */
 
void cgrdinit (
FILE *fp,
int  *ok
)
{
   int i, dummy,
       nres, natom, ncons;
 
   *ok = 1;
 
   rewind(fp);
 
   fscanf(fp,"%d %d %d",&nres,&natom,&ncons);
i = values.nres;
   if(values.nres != nres)
   {
      fprintf(out,"Warning from cgrdinit(): Current number of residues \
(%d) does not match\nnumber in conformations file (%d)\n",values.nres,nres);
      *ok = 0;
   }
 
   if(values.natoms != natom)
   {
      fprintf(out,"Warning from cgrdinit(): Current number of atom \
(%d) does not match\nnumber in conformations file (%d)\n",values.natoms,natom);
      *ok = 0;
   }
 
   /* Skip the constructed atoms list */
   for(i=0;i<ncons;i++)
      fscanf(fp,"%5d ",&dummy);
}
 

#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Read the constructed atom list from a CG file*/
 
void cgrdconsp (
FILE *fp,
int  *consp,
int  *ncons
)
{
   int nres,natom,i;
 
   rewind(fp);
 
   fscanf(fp,"%d %d %d",&nres,&natom,ncons);
 
   /* Read the constructed atoms list */
   for(i=0;i<*ncons;i++)
      fscanf(fp,"%5d ",&consp[i]);
}
 

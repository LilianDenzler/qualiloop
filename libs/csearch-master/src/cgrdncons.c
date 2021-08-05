#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Read number of constructed atoms*/
 
void cgrdncons (
FILE *fp,
int  *ncons
)
{
   rewind(fp);
   fscanf(fp,"%*d %*d %d",ncons);
}
 

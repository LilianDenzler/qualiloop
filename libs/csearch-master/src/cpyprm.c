#include "ProtoTypes.h"
#include "CongenProto.h"
 
void cpyprm(
float *parm,
short *itc,
short *iac,
float *atoma,
int *natom
)
{
   int i;
 
   if(*natom != 0)
      for(i=0; i< *natom; i++)
         atoma[i] = parm[itc[iac[i]-1]-1];
         /* FORTRAN: ATOMA(I)=PARM(ITC(IAC(I))) */
}
 

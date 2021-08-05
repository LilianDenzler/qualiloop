#include "ProtoTypes.h"
#include "CongenProto.h"
 
/************************************************************************
ACRM 20.06.91
SeaRches for WorD BounDed, i.e. searches the array A between LOW
and HI for WD. If successful, returns the FORTRAN index of WD in the
array. Otherwise, 0 is returned.
Uses pointers for low & high to be compatible with old FORTRAN version
************************************************************************/
 
int srwdbd(
char a[][4],
int *low,
int *high,
char wd[]
)
{
int i;
 
   if(*low)
   {
      for(i = *low - 1; i < *high; i++)
         if(!strncmp(a[i],wd,4))
            return(i+1);
   }
   return(0);
}
 

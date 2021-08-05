/* Sets up the triangular matrix parm_no, used for finding parameters */
#include "ProtoTypes.h"
#include "CongenProto.h"
 
void setup_parmno(
short parm_no[100][100]
)
{
int k=0,i,j;
 
   for(i=1;i<=100;i++)
   {
      for(j=1;j<=i;j++)
      {
         k++;
         parm_no[i-1][j-1] = k;
         parm_no[j-1][i-1] = k;
      }
   }
}

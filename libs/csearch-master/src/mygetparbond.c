#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* My C version of getparbond() */
 
float mygetparbond(
int *ib,
int *jb,
short int parm_no[100][100],
short int *iac,
int   *kcb,
int   *ncb,
char  type[maxat][4],
float *cbb
)
{
char temp1[8], temp2[8];
int i, j, icbb, code;
 
   if(*ib==0 || *jb==0)
      prdie("Error in mygetparbond(): Atom index is zero\n");
   i = iac[*ib-1];
   j = iac[*jb-1];
   code = (int)parm_no[j-1][i-1];
   icbb = bin_search(code,kcb,*ncb);
   if(icbb == 0)
   {
      strncpy(temp1,type[*ib-1],4);
      temp1[4] = '\0';
      strncpy(temp2,type[*jb-1],4);
      temp2[4] = '\0';
      fprintf(out,"Error in mygetparbond(): Unable to find parameters for \
bond between %4s and %4s\n",temp1,temp2);
      die();
   }
   return(cbb[icbb-1]);
}
 

#include "ProtoTypes.h"
#include "CongenProto.h"

#define WDMAX 80
 
/* Parses the restart string, and loads the degrees of freedom with
   the values found. 
   24.05.94 Corrected sscanf() call    By: ACRM
*/
 
void  setrestitr(
void
)
{
   struct dof *p;
   char wd[WDMAX];
   char *cp;
 
   cp = cg.restart_st;
   cp[cg.restart_stlen] = '\0';


   for (p=dof_head; p; p = p->next)
   {
      mynextwd(&cp,wd);
      p->restart_iter = 0;
      if (strlen(wd) > 0)
      {
/* ACRM 24.05.94
//        sscanf(wd,"%d",p->restart_iter);
*/
          sscanf(wd,"%d",&(p->restart_iter));
      }
   }
}
 

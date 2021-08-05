#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Set up an array of constructed atoms for all the degrees of freedom
*   from dof_head up to last_dof. If last_dof is null, then all degrees
*   of freedom are used.
*/
 
void setup_consp(
int **conspp,
int *nconspp,
struct dof *last_dof
)
{
   int *ip;
   struct dof *p;
   struct atom **atpp,*atp;
 
   *nconspp = 0;
   for (p=dof_head; p && (p->prev!=last_dof || last_dof == null); p = p->next)
      if (p->atoms)
         for (atpp=p->atoms; *atpp; atpp++)
            (*nconspp)++;
   ip = *conspp = (int *) alloc(*nconspp * sizeof(int));
   for (p=dof_head; p && (p->prev!=last_dof || last_dof == null); p = p->next)
      if (p->atoms)
         for (atpp=p->atoms; (atp = *atpp); atpp++)
            *ip++ = atp->atomno;
 
   /* ACRM changed to using sort_i() */
   sort_i((int *)*conspp,*nconspp);
}
 

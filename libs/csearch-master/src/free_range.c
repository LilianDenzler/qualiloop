#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Free the range list, rp. */
 
void free_range(
struct range_list *rp
)
{
   free(rp->low);
   free(rp->high);
   free(rp);
}
 

#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Calls print_atom1, used by FORPRTATM().
   Always sets the input string to a null string first. */
 
void c_print_atom(
char string[],
int  atnum
)
{
   string[0] = '\0';
   print_atom1(string,atnum);
}
 

#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Print the atoms pointed to by the array pointed to by atpp. */
 
void print_atoms(
struct atom **atpp
)
{
   fprintf(out,"Atoms for construction:\n");
   fprintf(out,"%s%s%s%s\n",
          "       Ante 0       ",
          "       Ante 1       ",
          "       Ante 2       ",
          "        Atom        ");
   while (*atpp) print_atom(*atpp++);
}
 

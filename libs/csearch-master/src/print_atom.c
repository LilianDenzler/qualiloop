#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Prints the antecedents and constructed atom for atp. */
 
void print_atom(
struct atom *atp
)
{
   char st[81];
 
   st[0] = '\0';
   print_atom1(st,atp->ante[0]);
   print_atom1(st,atp->ante[1]);
   print_atom1(st,atp->ante[2]);
   print_atom1(st,atp->atomno);
   fprintf(out,"%s\n",st);
}
 

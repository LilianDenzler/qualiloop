#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* This routine sets up for calculation fo CTER OT2 coords. */
 
void FixCterO(
PDB *pdb
)
{
char  temp[8];
 
PDB   *p,
      *ca_p = NULL,
      *c_p  = NULL,
      *o_p  = NULL;
 
   /* Step through the linked list */
   for(p=pdb;p;NEXT(p))
   {
      /* Record pointers to CA, C and O for each residue until we
         reach a CTER.      */
      if(!strncmp(p->resnam,"CTER",4) &&
         !strncmp(p->atnam,"OT2",3)   &&
         p->x == 9999.0)
         ccterxyz(p,ca_p,c_p,o_p);
 
      strncpy(temp,p->atnam,4);
      abmpad(temp,4);
      if(!strncmp(temp,"CA  ",4))
         ca_p = p;
      else if(!strncmp(temp,"C   ",4))
         c_p = p;
      else if(!strncmp(temp,"O   ",4) || !strncmp(temp,"OT1 ",4))
         o_p = p;
   }
}
 

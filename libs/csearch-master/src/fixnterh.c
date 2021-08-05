#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* This routine shuffles the pdb linked list as hadd() will have added
   Nter H's into the first true residue rather than into NTER */
 
void FixNterH(
PDB *pdb
)
{
   PDB   *p,
         *ht1_p,
         *ht2_p,
         *last_p;
   char  insert[8];
   int   found = 0,
         resnum;
 
   for(p=pdb;p;NEXT(p))
   {
      /* Look for NTER. When found, make a note of the atom pointers
         for the two H's and set the found flag.
      */
      if(!strncmp(p->resnam,"NTER",4))
      {
         last_p = p;
 
         if(!strncmp(p->atnam,"HT1",3))
         {
            found++;
            ht1_p = p;
         }
         else if(!strncmp(p->atnam,"HT2",3))
         {
            found++;
            ht2_p = p;
         }
      }
      else if(found==2 && strncmp(p->resnam,"NTER",4))
      {
         /* We're in the residue after an NTER and have found both H's.
            Step through this residue, finding HT1 and HT2, putting their
            coords into the NTER residue and deleting these atoms
         */
         for(resnum = p->resnum, strcpy(insert,p->insert);
             p && p->resnum==resnum && !strncmp(insert,p->insert,1);
             NEXT(p))
         {
            if(!strncmp(p->atnam,"HT1",3))
            {  /* Copy in the coords */
               ht1_p->x = p->x;
               ht1_p->y = p->y;
               ht1_p->z = p->z;
               /* Unlink this atom from the list */
               last_p->next = p->next;
               free(p);
               p = last_p;
            }
            else if(!strncmp(p->atnam,"HT2",3))
            {  /* Copy in the coords */
               ht2_p->x = p->x;
               ht2_p->y = p->y;
               ht2_p->z = p->z;
               /* Unlink this atom from the list */
               last_p->next = p->next;
               free(p);
               p = last_p;
            }
            else
            {
               last_p = p;
            }
         }
         found = 0;
      }
   }
}
 

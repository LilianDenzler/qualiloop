#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* AddH() merges a list of hydrogens for this atom into the main pdb
   structure list */
 
/* -- FUNCTION -- */
int AddH(
PDB *hlist,
PDB **position,
int igtype_m)
{
/* -- DEFINITIONS -- */
   PDB *p,*q,*r,*s;
   int atomcount=0 ;
   int k;
   int column ;
 
/* -- DECLARATIONS -- */
   extern char ghname[8][8] ;
 
/* -- CODE -- */
   q=hlist;
 
   /* Step through each atom in position list until we find the
      one corresponding to this PGP
   */
   for(k=1;k<16;k++)
   {
 
      p = position[k] ;
 
      if ( !p )
         column = 0 ;
      else if ( ( igtype_m==3 ) || ( igtype_m==5 ) )
/* For PGP types 3 & 5, look for atom in column 3 */
         column = 3 ;
      else if ( ( igtype_m==1 ) || ( igtype_m==2 ) || ( igtype_m==4 ) )
/* For PGP types 1,2,4 look in column 2 */
         column = 2 ;
      else
         column = 0 ;
 
/* Test the atom */
      if ( column && !(strncmp(p->atnam,ghname[column],4) ) )
      {
 
         /* Copy the atoms from hlist into the PDB list */
         s=p;
         r=p->next;           /* Store the pointer to the next record */
         ALLOCNEXT(p,PDB);    /* Insert a record in the main list */
         p->next=r;           /* Update its pointer */
         strcpy(p->junk,s->junk);   /* Copy the info into this record */
         p->atnum = ++atomcount;
         strcpy(p->atnam,q->atnam);
         strcpy(p->resnam,s->resnam);
         strcpy(p->chain,s->chain);
         p->resnum=s->resnum;
         strcpy(p->insert,s->insert);
         p->x=q->x;
         p->y=q->y;
         p->z=q->z;
         p->occ=1.0;
         p->bval=20.0;
 
         if((igtype_m==2)||(igtype_m==3))
         {
            /* For these types there are 2 or three H's to add, so
               add another one
            */
            NEXT(q);
            s=p;
            r=p->next;
            ALLOCNEXT(p,PDB);
            p->next=r;
            strcpy(p->junk,s->junk);
            p->atnum = ++atomcount;
            strcpy(p->atnam,q->atnam);
            strcpy(p->resnam,s->resnam);
            strcpy(p->chain,s->chain);
            p->resnum=s->resnum;
            strcpy(p->insert,s->insert);
            p->x=q->x;
            p->y=q->y;
            p->z=q->z;
            p->occ=1.0;
            p->bval=20.0;
         }
         if(igtype_m==3)
         {
            /* For this type there are 3 H's to add, so
               add another one
            */
            NEXT(q);
            s=p;
            r=p->next;
            ALLOCNEXT(p,PDB);
            p->next=r;
            strcpy(p->junk,s->junk);
            p->atnum = ++atomcount;
            strcpy(p->atnam,q->atnam);
            strcpy(p->resnam,s->resnam);
            strcpy(p->chain,s->chain);
            p->resnum=s->resnum;
            strcpy(p->insert,s->insert);
            p->x=q->x;
            p->y=q->y;
            p->z=q->z;
            p->occ=1.0;
            p->bval=20.0;
         }
      return(0);
      }   /* End of matches */
   }  /* End of main list */
   return(0);
}
 

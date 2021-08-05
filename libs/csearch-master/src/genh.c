#include "OMLConstants.h"
#include "ProtoTypes.h"
#include "CongenProto.h"

extern int info_level;

/* Returns the number of hydrogens created, or 0, if an error */
 
int GenH(
PDB *pdb,
unsigned short *err_flag,
char grname[8]
)
{
   unsigned short  firstres;
   char  *bl       = "    ",
         *co       = "CO  ",
         *c        = "C   ";
   int   k, n, m, j, ittot;
   PDB   *p,*q;
   PDB   *position[MAXATINRES],*hlist;
 
   for(j=0;j<MAXATINRES;j++)
      position[j]=NULL;
 
   hfac    = 180.0/PI;
   fac     = 0.5 * sqrt(3.0);
 
   firstres = TRUE;
 
   /* Main loop */
   ih=1; it1=0; it2=0; it3=0; it4=0; it5=0;
 
   /* For the first residue, we don't have info for the previous C
      so set gnat[1] (the atom list for this residue) to a blank */
   strcpy(gnat[1],bl);
 
   position[1]=NULL;
 
   /* For each atom in the PDB file */
   for(p=pdb;p;)
   {
      /* Don't do anything with NTER residues except reset firstres */
      if(!strncmp(p->resnam,"NTER",4))
      {
         firstres = TRUE;
         NEXT(p);
         continue;
      }
 
      k=1;  /* This is one as this position is reserved for the CO of
               the previous residue */
 
      /* Copy this residue into our global work arrays.
         grname         is the residue name
         gx[],gy[],gz[] are the coordinates of the atoms
         gnat[]         are the atom names      */
      do
      {
         k++;
         position[k] = p;
 
         no=p->resnum;
         strcpy(grname,p->resnam);
#ifdef DEBUG
         printf("Group name is: %s\n",grname);
#endif
         strcpy(gnat[k],p->atnam);
#ifdef DEBUG
         printf("Atom %d name: %s\n",k,gnat[k]);
#endif
         gx[k] = p->x;
         gy[k] = p->y;
         gz[k] = p->z;
         NEXT(p);
         if(!p) break;
      } while(p && p->resnum == no);
 
      /* kmax is used to store the number of atoms in this residue */
      kmax = k;
 
      /* Go through each of the PGP's until we find this residue type */
      for(n=1;n<=npgp;n++)
      {
#ifdef DEBUG
         printf("gres(%d)=%s\n",n,gres[n]);
#endif
         if(strncmp(grname,gres[n],4)) continue;
#ifdef DEBUG
         printf("Entry found for this res. in PGP\n");
#endif
 
         /* Having found the residue type, copy the associated PGP atom
            list into the ghname[] work array */
         for(m=1;m<=6;m++) strcpy(ghname[m],gatom[n][m]);
 
         /* Now generate the hydrogen associated with this PGP */
         if((hlist = makeh(igtype[n],gr[n],alpha[n],beta[n],firstres))!=NULL)
            AddH(hlist,position,igtype[n]);
            /* and add it into the list, updating p to point to the new
               end of this residue */
         if(*err_flag) return(0);
      }
 
      /* If this is the first residue then handle it as NTER */
      if(firstres)
      {
#ifdef DEBUG
         printf("Rechecking for NTER...\n");
#endif
         for(n=1;n<=npgp;n++)
         {
#ifdef DEBUG
            printf("gres(%d)=%s\n",n,gres[n]);
#endif
            if(strncmp("NTER",gres[n],4)) continue;
#ifdef DEBUG
            printf("Entry found for NTER in PGP.\n");
#endif
 
            /* Having found the residue type, copy the associated PGP atom
               list into the ghname[] work array */
            for(m=1;m<=6;m++) strcpy(ghname[m],gatom[n][m]);
 
            /* Now generate the hydrogen associated with this PGP */
            if((hlist = makeh(igtype[n],gr[n],alpha[n],beta[n],firstres))!=NULL)
               AddH(hlist,position,igtype[n]);
               /* and add it into the list, updating p to point to the new
                  end of this residue */
            if(*err_flag) return(0);
         }
      }
 
      /* Next amino acid
         Set up pointer for carbonyl of previous residue unless this
         residue is a CTER when we set firstres to TRUE      */
 
      if(strncmp(grname,"CTER",4))
      {
         firstres = FALSE;
         for(j=kmax-1;j>0;j--) if(!strncmp(gnat[j],c,4))break;
         if(j==0)
         {
            printf("\nError==> genh() found no carbonyl carbon\
 in residue %d\n\n",p->resnum);
            *err_flag=1;
            return(0);
         }
         strcpy(gnat[1],co);
         gx[1]=gx[j];
         gy[1]=gy[j];
         gz[1]=gz[j];
         q=position[j];
 
         for(k=0;k<MAXATINRES;position[k++]=NULL);
         position[1]=q;
      }
      else
      {
         firstres = TRUE;
      }
 
   }  /* Go back to the next atom/residue */
 
   ittot=it1+it2+it3+it4+it5;
 
   if(info_level>0)
      printf("Total number of protons generated  =%5d\n\
             Type 1 C-H's          =%5d\n\
             Type 2 C-H2's         =%5d\n\
             Type 3 C-H3's         =%5d\n\
             Type 4 sp2 C-H's,>N-H =%5d\n\
             Type 5 O-H's =N-H's   =%5d\n",
             ittot,it1,it2,it3,it4,it5);
 
   return(ittot);
}
 

/* #define SCREEN_INFO */  /* Causes error messages to be output by screen()  */
/*************************************************************************

   Program:    
   File:       HAddPDB.c
   
   Version:    V2.3
   Date:       27.07.93
   Function:   
   
   Copyright:  SciTech Software 1991
   Author:     Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0372) 275775
   EMail:      UUCP: cbmuk!cbmuka!scitec!amartin
               JANET: andrew@uk.ac.ox.biop
               
   Original version written while at:
               Laboratory of Mathematical Biology
               National Institute for Medical Research
               The Ridgeway
               Mill Hill
               London
               NW7 1AA
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission from
   the author, although it may be given away free with commercial products,
   providing it is made clear that this program is free and that the source
   code is provided with the program.

**************************************************************************

   Description:
   ============
   Routine to add hydrogens to a protein linked list of type PDB.
   The routine allocates space for the new atoms and inserts them
   into the list at the appropriate positions within the residues.

**************************************************************************

   Usage:
   ======
   hadd(fp,pdb)
   Input:         FILE   *fp        File containing proton generation
                                    parameters.
   Input/Output:  PDB    *pdb       Linked list of protein structure.
   Returns:       int               Number of hydrogens added.
   Externs:       int    info_level Message level.

**************************************************************************

   Revision History:
   =================
   V2.0  16.05.90 AddH is changed to insert each set of atoms for each 
                  PGP, on the fly, rather than building a complete list 
                  of hydrogens and then merging the two lists. This allows
                  us to get round the problem of missing atoms, since 
                  there will be no merging error.

   V2.1  24.05.90 Returns the number fo hydrogens added. Also fixes bug 
                  relating to number of type 2 and type 3 H's added.
                  Doesn't work under UNIX!

   V2.2  15.07.91 A few bits of tidying up:
                  >  Now uses macros.h rather than defining macros itself. 
                  >  Arrays now changed so should fix alignment problems 
                     under UNIX. 
                  >  Now correctly checks return from forscanf() and no 
                     longer reads characters to check EOF itself.
                  >  Improves treatment of NTER residues where it now 
                     generates H's on the first true residue. Residues 
                     labelled NTER will be ignored. Calling FixNterH() 
                     will move the H coords into the NTER residue if 
                     required. 
                  >  Fixes reading of PGP files with blank lines 
                  >  Bug fix for type 3's
                  Currently untested under UNIX.

   V2.3  27.07.93 Changed to use fsscanf() and I/O precision is double

*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "bioplib/fsscanf.h"
#include "bioplib/macros.h"

#ifdef SCREEN_INFO
#include "bioplib/WindIO.h"
#endif

#define MAXTYPE 200

#define CLEAR_PDB(p) strcpy(p->junk,"      "); \
                     p->atnum=0; \
                     strcpy(p->atnam,"    "); \
                     strcpy(p->resnam,"    "); \
                     p->resnum=0; \
                     strcpy(p->insert," "); \
                     strcpy(p->chain," "); \
                     p->x = 0.0; p->y = 0.0; p->z = 0.0; \
                     p->occ = 0.0; p->bval = 0.0; \
                     p->next = NULL

/************************************************************************/
/* Global Arrays, etc.
*/
static char    gres[MAXTYPE][8],
               gatom[MAXTYPE][8][8],
               grname[8],
               ghname[8][8],
               gnat[16][8];
static int     igtype[MAXTYPE],
               npgp,
               no,
               kmax,
               it1, it2, it3, it4, it5, ih;
static REAL    gr[MAXTYPE],
               alpha[MAXTYPE],
               beta[MAXTYPE],
               gx[16], gy[16], gz[16],
               hfac, fac;

/* This one is externally visible                                       */
struct
{
   int   Total,      /* Total hydrogens                                 */
         T1,         /* Type 1 C-H's                                    */
         T2,         /* Type 2 C-H2's                                   */
         T3,         /* Type 3 C-H3's                                   */
         T4,         /* Type 4 sp2 C-H's,>N-H                           */
         T5;         /* Type 5 O-H's =N-H's                             */
}  HaddInfo;


/************************************************************************/
/* Prototypes
*/
/* #include "HAddPDB.p" */

int HAddPDB(FILE *fp, PDB  *pdb);
int ReadPGP(FILE *fp);
int GenH(PDB *pdb, unsigned short *err_flag);
PDB *makeh(int    igtype_m,
           REAL   gr_m,
           REAL   alpha_m,
           REAL   beta_m,
           BOOL   firstres);
BOOL AddH(PDB *hlist, PDB **position, int igtype_m);

/************************************************************************/
/*>int HAddPDB(FILE *fp, PDB  *pdb)
   -----------------------------
   This is the main entry point. Returns number of hydrogens added, 0 if
   error
   04.01.94 Changed check on return=NULL to 0
*/
int HAddPDB(FILE *fp, PDB  *pdb)
{
   PDB *p;
   int atomcount=1;
   int nhydrogens;
   unsigned short err_flag=0;

   /* Read the parameter file                                           */
   ReadPGP(fp);

   /* Generate the hydrogens                                            */
   if((nhydrogens=GenH(pdb,&err_flag))==0)
   {
      return(0);
   }
   
   /* Renumber atoms in PDB linked list                                 */
   for(p=pdb;p;NEXT(p)) p->atnum=atomcount++;
   return(nhydrogens);
}



/************************************************************************/
/*>int ReadPGP(FILE *fp)
   ---------------------
   Read a proton generation parameter file
*/
int ReadPGP(FILE *fp)
{
   char  buffer[160];
   int   n=0;
   
   while(fgets(buffer,159,fp))
   {
      fsscanf(buffer,
              "%4s%4s%1x%4s%1x%4s%1x%4s%1x%4s%1x%4s%1x%1d%10lf%10lf%10lf",
              gres[++n],              
              gatom[n][1],
              gatom[n][2],
              gatom[n][3],
              gatom[n][4],
              gatom[n][5],
              gatom[n][6],
              &(igtype[n]),
              &(gr[n]),
              &(alpha[n]),
              &(beta[n]));

#ifdef DEBUG
      printf("%4s %4s %4s %4s %4s %4s %4s %1d %8.3lf %8.3lf %8.3lf\n",
             gres[n],gatom[n][1],gatom[n][2],gatom[n][3],gatom[n][4],
             gatom[n][5],gatom[n][6],igtype[n],gr[n],alpha[n],beta[n]);
#endif

      if(igtype[n] != 0)
      {
         alpha[n] *= (PI/180.0);    
         beta[n]  *= (PI/180.0);
      }
      
   }  /* End of file                                                    */

   npgp = n;

   return(npgp);
}


/************************************************************************/
/*>int GenH(PDB *pdb, unsigned short *err_flag)
   --------------------------------------------
   Returns the number of hydrogens created, or 0, if an error
*/
int GenH(PDB *pdb, unsigned short *err_flag)
{
   BOOL  firstres;
   char  *bl       = "    ",
         *co       = "CO  ",
         *c        = "C   ";
   int   k, n, m, j, ittot;
   PDB   *p,*q;
   PDB   *position[16],*hlist;
#ifdef SCREEN_INFO
   char  buffer[160];
#endif
    
   for(j=0;j<16;position[j++]=NULL);

   hfac    = 180.0/PI,
   fac     = 0.5 * sqrt(3.0);
    
   firstres = TRUE;
    
   /* Main loop                                                         */
   ih=1; it1=0; it2=0; it3=0; it4=0; it5=0;

   /* For the first residue, we don't have info for the previous C
      so set gnat[1] (the atom list for this residue) to a blank        */
   strcpy(gnat[1],bl);
   
   position[1]=NULL;

   /* For each atom in the PDB file                                     */
   for(p=pdb;p;)
   {
      /* Don't do anything with NTER residues except reset firstres     */
      if(!strncmp(p->resnam,"NTER",4))
      {
         firstres = TRUE;
         NEXT(p);
         continue;
      }

      k=1;  /* This is one as this position is reserved for the CO of
               the previous residue 
            */

      /* Copy this residue into our global work arrays.
         grname         is the residue name
         gx[],gy[],gz[] are the coordinates of the atoms
         gnat[]         are the atom names
      */
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
      } while(p,p->resnum == no);

      /* kmax is used to store the number of atoms in this residue      */
      kmax = k;

      /* Go through each of the PGP's until we find this residue type   */
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
            list into the ghname[] work array
         */
         for(m=1;m<=6;m++) strcpy(ghname[m],gatom[n][m]);
         
         /* Now generate the hydrogen associated with this PGP          */
         if((hlist = makeh(igtype[n],gr[n],alpha[n],beta[n],firstres))
            !=NULL)
         {
            /* And add it into the list, updating p to point to the new
               end of this residue
            */
            if(!AddH(hlist,position,igtype[n])) return(0);
         }
         if(*err_flag) return(0);
      }
      
      /* If this is the first residue then handle it as NTER            */
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
               list into the ghname[] work array                        */
            for(m=1;m<=6;m++) strcpy(ghname[m],gatom[n][m]);
         
            /* Now generate the hydrogen associated with this PGP       */
            if((hlist = makeh(igtype[n],gr[n],alpha[n],beta[n],firstres))
               !=NULL)
            {
               /* And add it into the list, updating p to point to the new
                  end of this residue 
               */
               if(!AddH(hlist,position,igtype[n])) return(0);
            }
            if(*err_flag) return(0);
         }
      }

      /* Next amino acid
         Set up pointer for carbonyl of previous residue unless this 
         residue is a CTER when we set firstres to TRUE
      */
      
      if(strncmp(grname,"CTER",4))
      {
         firstres = FALSE;
         for(j=kmax-1;j>0;j--) if(!strncmp(gnat[j],c,4))break;
         if(j==0)
         {
#ifdef SCREEN_INFO
            sprintf(buffer,"\nError==> genh() found no carbonyl carbon \
in residue %d\n\n",p->resnum);
            screen(buffer);
#endif
            *err_flag=1;
            return(0);
         }
         strcpy(gnat[1],co);
         gx[1]=gx[j];
         gy[1]=gy[j];
         gz[1]=gz[j];
         q=position[j];

         for(k=0;k<16;position[k++]=NULL);
         position[1]=q;
      }
      else
      {
         firstres = TRUE;
      }
      
   }  /* Go back to the next atom/residue                               */
   
   ittot=it1+it2+it3+it4+it5;
   
   HaddInfo.Total = ittot;
   HaddInfo.T1    = it1;
   HaddInfo.T2    = it2;
   HaddInfo.T3    = it3;
   HaddInfo.T4    = it4;
   HaddInfo.T5    = it5;

   return(ittot);
}

/************************************************************************/
/*>PDB *makeh(int igtype_m, REAL gr_m, REAL alpha_m, REAL beta_m, 
              BOOL firstres)
   --------------------------------------------------------------
   Generate a hydrogen coordinate. Returns NULL if fails to allocate
   memory or this is the Nter N where we don't require a planar H or
   all atoms have been done.
*/
PDB *makeh(int    igtype_m,
           REAL   gr_m,
           REAL   alpha_m,
           REAL   beta_m,
           BOOL   firstres)
{ 
   char  *nt = "NT  ",
         *n  = "N   ";
   REAL  x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,
         x21,y21,z21,r21,
         x21p,y21p,z21p,r21p,
         x23,y23,z23,r23,
         xp23,yp23,zp23,rp23,
         x24,y24,z24,r24,
         x32,y32,z32,
         xv25,yv25,zv25,rv25,
         cosa,sina,
         cosb,sinb,
         cosax,cosay,cosaz,
         xa,ya,za,xb,yb,zb,
         xab,yab,zab,rab,
         xmin,ymin,zmin,
         xapb,yapb,zapb,rapb,
         xplus,yplus,zplus,
         xp,yp,zp,
         xs,ys,zs,
         xh,yh,zh,
         xv,yv,zv,
         scalpr;
   int   kount = 0,
         num_ant,k,jj;
   unsigned short nt_point,ok;
   PDB   *hlist,*hlist_start;

   INIT(hlist_start, PDB);
   if(hlist_start == NULL) return(NULL);
   
   CLEAR_PDB(hlist_start);
   hlist = hlist_start;

   nt_point = 0;
   if(igtype_m==1) num_ant=4; else num_ant=3;

   /* Don't add a planar H to the Nter N.                               */
   if(firstres && igtype_m==4 && !strncmp(ghname[2],n,4)) return(NULL);

   /* Work through the atoms in this residue (gnat[]) and compare them 
      with the first 4 atoms in the PGP atom list in ghname[], storing 
      the associated coordinates
   */
   for(k=1;k<=kmax;k++)
   {
      if(nt_point) 
      {
         strcpy(gnat[nt_point],"N   ");
         nt_point=0;
      }

      if(!strncmp(ghname[1],gnat[k],4))
      {
         if(!strncmp(gnat[k],"NT  ",4)) nt_point=k;
         kount++;
         x1=gx[k];
         y1=gy[k];
         z1=gz[k];
      }
      else if(!strncmp(ghname[2],gnat[k],4))
      {
         if(!strncmp(gnat[k],"NT  ",4)) nt_point=k;
         kount++;
         x2=gx[k];
         y2=gy[k];
         z2=gz[k];
      }
      else if(!strncmp(ghname[3],gnat[k],4))
      {
         if(!strncmp(gnat[k],"NT  ",4)) nt_point=k;
         kount++;
         x3=gx[k];
         y3=gy[k];
         z3=gz[k];
      }
      else if(!strncmp(ghname[4],gnat[k],4))
      {
         if(!strncmp(gnat[k],"NT  ",4)) nt_point=k;
         kount++;
         x4=gx[k];
         y4=gy[k];
         z4=gz[k];
      }
   }  /* End of k loop around this residue                              */

   /* Check we found all the atoms we need for this PGP                 */
   if(kount != num_ant)
   {
      ok = FALSE;
      if(firstres)
      {
         /* Check it's not the missing N in the first residue           */
         for(jj=1;jj<=4;jj++) if(!strncmp(ghname[jj],n,4)) ok = TRUE;
      }
      else
      {
         /* Check it's not just the NT                                  */
         for(jj=1;jj<=4;jj++) if(!strncmp(ghname[jj],nt,4)) ok = TRUE;
      }
#ifdef SCREEN_INFO
      if(!ok)
      {
         char buffer[160];
         
         sprintf(buffer,"Error==> makeh() unable to find all atoms \
required by PGP parameter for %3s %5d\n",grname,no);
         screen(buffer);
         screen("Atoms required by PGP\n");
         screen("GHNAME: ");
         for(jj=1;jj<=4;jj++)
         {
            sprintf(buffer," %4s",ghname[jj]);
            screen(buffer);
         }
         screen("\n");
         screen("Atoms in current residue\n");
         screen("GNAT  : ");
         for(jj=1;jj<=kmax;jj++)
         {
            sprintf(buffer," %4s",gnat[jj]);
            screen(buffer);
         }
         screen("\n");
      }
#endif
      return(NULL);
   }
    
   x21=x2-x1;
   y21=y2-y1;
   z21=z2-z1;
   r21=sqrt(x21*x21 + y21*y21 + z21*z21);

   x23=x2-x3;
   y23=y2-y3;
   z23=z2-z3;
   r23=sqrt(x23*x23 + y23*y23 + z23*z23);
    
/* IGTYPE 1: Generation of 1 tetrahedral H                              */
    
   if(igtype_m == 1)
   {
      x24=x2-x4;
      y24=y2-y4;
      z24=z2-z4;
      r24=sqrt(x24*x24 + y24*y24 + z24*z24);
      xv25=x21/r21+x24/r24+x23/r23;
      yv25=y21/r21+y24/r24+y23/r23;
      zv25=z21/r21+z24/r24+z23/r23;
      rv25=sqrt(xv25*xv25 + yv25*yv25 + zv25*zv25);
      x5=x2+gr_m*xv25/rv25;
      y5=y2+gr_m*yv25/rv25;
      z5=z2+gr_m*zv25/rv25;

      hlist->resnum=no;
      strcpy(hlist->atnam,ghname[5]);
      hlist->x=x5;
      hlist->y=y5;
      hlist->z=z5;
      ALLOCNEXT(hlist,PDB);
      if(hlist == NULL)
      {
         FREELIST(hlist_start, PDB);
         return(NULL);
      }
      CLEAR_PDB(hlist);

#ifdef DEBUG
      printf("makeh() Type 1 allocated hlist = %d\n",(int)hlist);
#endif

      it1++;
   }  /* End of IGTYPE 1                                                */
   else
   {
      cosa=cos(alpha_m);
      sina=sin(alpha_m);
      switch(igtype_m)
      {
/* IGTYPE 2: Generation of 2 tetrahedral H's                            */
case 2:  xa=x21/r21;
         ya=y21/r21;
         za=z21/r21;
         xb=x23/r23;
         yb=y23/r23;
         zb=z23/r23;
         xab=xa-xb;
         yab=ya-yb;
         zab=za-zb;
         rab=sqrt(xab*xab+yab*yab+zab*zab);
         xmin=xab/rab;
         ymin=yab/rab;
         zmin=zab/rab;
         xapb=xa+xb;
         yapb=ya+yb;
         zapb=za+zb;
         rapb=sqrt(xapb*xapb+yapb*yapb+zapb*zapb);
         xplus=xapb/rapb;
         yplus=yapb/rapb;
         zplus=zapb/rapb;
         xs=yplus*zmin-zplus*ymin;
         ys=zplus*xmin-xplus*zmin;
         zs=xplus*ymin-yplus*xmin;
         x4=x2+gr_m*(cosa*xplus+sina*xs);
         y4=y2+gr_m*(cosa*yplus+sina*ys);
         z4=z2+gr_m*(cosa*zplus+sina*zs);
         x5=x2+gr_m*(cosa*xplus-sina*xs);
         y5=y2+gr_m*(cosa*yplus-sina*ys);
         z5=z2+gr_m*(cosa*zplus-sina*zs);

         hlist->resnum=no;
         strcpy(hlist->atnam,ghname[4]);
         hlist->x=x4;
         hlist->y=y4;
         hlist->z=z4;
         ALLOCNEXT(hlist,PDB);
         if(hlist == NULL)
         {
            FREELIST(hlist_start, PDB);
            return(NULL);
         }
         CLEAR_PDB(hlist);

#ifdef DEBUG
         printf("makeh() Type 2a allocated hlist = %d\n",(int)hlist);
#endif

         hlist->resnum=no;
         strcpy(hlist->atnam,ghname[5]);
         hlist->x=x5;
         hlist->y=y5;
         hlist->z=z5;
         ALLOCNEXT(hlist,PDB);
         if(hlist == NULL)
         {
            FREELIST(hlist_start, PDB);
            return(NULL);
         }
         CLEAR_PDB(hlist);

#ifdef DEBUG
         printf("makeh() Type 2b allocated hlist = %d\n",(int)hlist);
#endif

         it2+=2;
         break;

/* Initialisation for both these cases is the same                      */
case 3:
case 5:  x32=x3-x2;
         y32=y3-y2;
         z32=z3-z2;
         xh=x32/r23;
         yh=y32/r23;
         zh=z32/r23;
         scalpr=(x21*x32+y21*y32+z21*z32)/r23;
         xp=scalpr*xh;
         yp=scalpr*yh;
         zp=scalpr*zh;
         x21p=x21-xp;
         y21p=y21-yp;
         z21p=z21-zp;
         r21p=sqrt(x21p*x21p+y21p*y21p+z21p*z21p);
         xv=x21p/r21p;
         yv=y21p/r21p;
         zv=z21p/r21p;
         xs=yh*zv-zh*yv;
         ys=zh*xv-xh*zv;
         zs=xh*yv-yh*xv;
         cosax=cosa*xh;
         cosay=cosa*yh;
         cosaz=cosa*zh;

/* IGTYPE 3: Generation of 3 tetrahedral H's                            */
         if(igtype_m==3)
         {
            x4=x3+gr_m*(cosax+sina*xv);
            y4=y3+gr_m*(cosay+sina*yv);
            z4=z3+gr_m*(cosaz+sina*zv);
            
            /* V2.2: Bug fix here: xy, ys, zs; not xs all the time!     */
            x5=x3+gr_m*(cosax+sina*(fac*xs-0.5*xv));
            y5=y3+gr_m*(cosay+sina*(fac*ys-0.5*yv));
            z5=z3+gr_m*(cosaz+sina*(fac*zs-0.5*zv));
            x6=x3+gr_m*(cosax+sina*(-fac*xs-0.5*xv));
            y6=y3+gr_m*(cosay+sina*(-fac*ys-0.5*yv));
            z6=z3+gr_m*(cosaz+sina*(-fac*zs-0.5*zv));

            hlist->resnum=no;
            strcpy(hlist->atnam,ghname[4]);
            hlist->x=x4;
            hlist->y=y4;
            hlist->z=z4;
            ALLOCNEXT(hlist,PDB);
            if(hlist == NULL)
            {
               FREELIST(hlist_start, PDB);
               return(NULL);
            }
            CLEAR_PDB(hlist);

#ifdef DEBUG
            printf("makeh() Type 3a allocated hlist = %d\n",(int)hlist);
#endif

            hlist->resnum=no;
            strcpy(hlist->atnam,ghname[5]);
            hlist->x=x5;
            hlist->y=y5;
            hlist->z=z5;
            ALLOCNEXT(hlist,PDB);
            if(hlist == NULL)
            {
               FREELIST(hlist_start, PDB);
               return(NULL);
            }
            CLEAR_PDB(hlist);

#ifdef DEBUG
            printf("makeh() Type 3b allocated hlist = %d\n",(int)hlist);
#endif

            hlist->resnum=no;
            strcpy(hlist->atnam,ghname[6]);
            hlist->x=x6;
            hlist->y=y6;
            hlist->z=z6;
            ALLOCNEXT(hlist,PDB);
            if(hlist == NULL)
            {
               FREELIST(hlist_start, PDB);
               return(NULL);
            }
            CLEAR_PDB(hlist);

#ifdef DEBUG
            printf("makeh() Type 3c allocated hlist = %d\n",(int)hlist);
#endif

            it3+=3;
         }
         else if(igtype_m==5)
         {
/* IGTYPE 5: Generation of 1 H where an angle is specified              */
            cosb=cos(beta_m);
            sinb=sin(beta_m);
            x4=x3+gr_m*(cosax+sina*(cosb*xv+sinb*xs));
            y4=y3+gr_m*(cosay+sina*(cosb*yv+sinb*ys));
            z4=z3+gr_m*(cosaz+sina*(cosb*zv+sinb*zs));

            hlist->resnum=no;
            strcpy(hlist->atnam,ghname[4]);
            hlist->x=x4;
            hlist->y=y4;
            hlist->z=z4;
            ALLOCNEXT(hlist,PDB);
            if(hlist == NULL)
            {
               FREELIST(hlist_start, PDB);
               return(NULL);
            }
            CLEAR_PDB(hlist);

#ifdef DEBUG
            printf("makeh() Type 5 allocated hlist = %d\n",(int)hlist);
#endif

            it5++;
         }
         break;

/* IGTYPE 4: Generation of 1 sp2 H                                      */
case 4:  x32=x3-x2;
         y32=y3-y2;
         z32=z3-z2;
         scalpr=(x21*x32+y21*y32+z21*z32)/r21;
         xh=x21/r21;
         yh=y21/r21;
         zh=z21/r21;
         xp=scalpr*xh;
         yp=scalpr*yh;
         zp=scalpr*zh;
         xp23=xp+x23;
         yp23=yp+y23;
         zp23=zp+z23;
         rp23=sqrt(xp23*xp23+yp23*yp23+zp23*zp23);
         xv=xp23/rp23;
         yv=yp23/rp23;
         zv=zp23/rp23;
         x4=x2+gr_m*(sina*xv-cosa*xh);
         y4=y2+gr_m*(sina*yv-cosa*yh);
         z4=z2+gr_m*(sina*zv-cosa*zh);

         hlist->resnum=no;
         strcpy(hlist->atnam,ghname[4]);
         hlist->x=x4;
         hlist->y=y4;
         hlist->z=z4;
         ALLOCNEXT(hlist,PDB);
         if(hlist == NULL)
         {
            FREELIST(hlist_start, PDB);
            return(NULL);
         }
         CLEAR_PDB(hlist);

#ifdef DEBUG
         printf("makeh() Type 4 allocated hlist = %d\n",(int)hlist);
#endif

         it4++;
      }  /* End of switch                                               */
   }  /* End of IGTYPE 1 else clause                                    */

   return(hlist_start);
}

/************************************************************************/
/*>BOOL AddH(PDB *hlist, PDB **position, int igtype_m)
   ---------------------------------------------------
   AddH() merges a list of hydrogens for this atom into the main pdb 
   structure list. Returns FALSE if the procedure failed.
*/
BOOL AddH(PDB *hlist, PDB **position, int igtype_m)
{
   PDB *p,*q,*r,*s;
   int atomcount=0,
       k;

   q=hlist;

   /* Step through each atom in position list until we find the
      one corresponding to this PGP
   */
   for(k=1;k<16;k++)
   {
      if(!position[k]) continue;
      p = position[k];

      /* For PGP types 3 & 5, look for atom in column 3                 */
      if((((igtype_m==3)||(igtype_m==5))
        &&(!(strncmp(p->atnam,ghname[3],4))))||
        /* For PGP types 1,2,4 look in column 2                         */
        (((igtype_m==1)||(igtype_m==2)||(igtype_m==4))
        &&(!strncmp(p->atnam,ghname[2],4))))
      {

         /* Copy the atoms from hlist into the PDB list                 */
         s=p;
         r=p->next;           /* Store the pointer to the next record   */
         ALLOCNEXT(p,PDB);    /* Insert a record in the main list       */
         if(p == NULL)
            return(FALSE);
         
         p->next=r;           /* Update its pointer                     */
         strcpy(p->junk,s->junk);   /* Copy the info into this record   */
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
            if(p==NULL)
               return(FALSE);
               
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
            if(p==NULL)
               return(FALSE);
               
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
      }   /* End of matches                                             */
   }  /* End of main list                                               */
   return(TRUE);
}

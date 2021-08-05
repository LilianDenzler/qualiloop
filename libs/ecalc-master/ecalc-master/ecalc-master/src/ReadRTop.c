/*************************************************************************

   Program:    ECalc
   File:       ReadRTop.c
   
   Version:    V1.5.2
   Date:       05.03.21
   Function:   Read a residue topology file into global arrays
   
   Copyright:  (c) UCL, Prof. Andrew C. R. Martin 1994-2021
   Author:     Prof. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

   Reads the parameter file into global arrays.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1   30.08.94 Original
   V0.2   22.09.94 Added code to get file from environment variable
   V1.0   30.09.94 First release version
   V1.1   12.10.94 Changes to energy.c
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Fixed for new version of GetWord()
   V1.5.1 07.01.21 Skipped
   V1.5.2 05.03.21 Skipped

*************************************************************************/
/* Includes
*/
#define READRTOP_MAIN

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/parse.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

#include "rtop.h"
#include "ecalc.h"

/************************************************************************/
/* Defines and macros
*/
#define KEY_RESIDUE  0
#define KEY_ATOM     1
#define KEY_BOND     2
#define KEY_ANGLE    3
#define KEY_TORSION  4
#define KEY_IMPROPER 5
#define KEY_DONOR    6
#define KEY_ACCEPTOR 7
#define KEY_BUILD    8
#define NKEYS        9
#define MAXSTRPARAM  10
#define MAXREALPARAM 1
#define MAXSTRLEN    80  /* Fixed from 16(!)  14.09.94                  */

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static BOOL StoreResidueRTop(char *res, REAL charge);
static BOOL StoreAtomRTop(char *atom, char *type, REAL charge, 
                          char *exclusions);
static BOOL StoreBondRTop(char *atom1, char *atom2);
static BOOL StoreAngleRTop(char *atom1, char *atom2, char *atom3);
static BOOL StoreTorsionRTop(char *atom1, char *atom2, 
                             char *atom3, char *atom4);
static BOOL StoreImproperRTop(char *atom1, char *atom2, 
                              char *atom3, char *atom4);
static BOOL StoreDonorRTop(char *atom1, char *atom2, 
                           char *atom3, char *atom4);
static BOOL StoreAcceptorRTop(char *atom);

/************************************************************************/
/*>BOOL ReadRTop(char *filename)
   -----------------------------
   Sets up the parser and parses the residue topology file dispatching
   routines to fill in the data.

   30.08.94 Original    By: ACRM
   14.09.94 Uses sscanf instead of atof()
   22.09.94 Changed to call OpenFile() rather than fopen()
*/
BOOL ReadRTop(char *filename)
{
   KeyWd keywords[NKEYS];
   char  *StrParam[MAXSTRPARAM],
         buffer[MAXBUFF];
   REAL  RealParam[MAXREALPARAM],
         temp;
   int   i;
   FILE  *fp;
   BOOL  NoVar;
   
   /* Open the file                                                     */
   if((fp=OpenFile(filename,DATADIR,"r",&NoVar)) == NULL)
   {
      sprintf(gError,"Unable to open file `%s'",filename);
      StoreError("ReadRTop()",gError);

      if(NoVar)
      {
         sprintf(gError,"Environment variable `%s' not set",DATADIR);
         StoreError("ReadRTop()",gError);
      }
      
      return(FALSE);
   }

   /* Allocate memory for return strings                                */
   for(i=0; i<MAXSTRPARAM; i++)
   {
      if((StrParam[i] = (char *)malloc(MAXSTRLEN * sizeof(char))) == NULL)
      {
         StoreError("ReadRTop()","No memory for parser");
         return(FALSE);
      }
   }
   
   /* Build keywords                                                    */
   MAKEKEY(keywords[KEY_RESIDUE],  "RESIDUE",  STRING, 2);
   MAKEKEY(keywords[KEY_ATOM],     "ATOM",     STRING, 4);
   MAKEKEY(keywords[KEY_BOND],     "BOND",     STRING, 2);
   MAKEKEY(keywords[KEY_ANGLE],    "ANGLE",    STRING, 3);
   MAKEKEY(keywords[KEY_TORSION],  "TORSION",  STRING, 4);
   MAKEKEY(keywords[KEY_IMPROPER], "IMPROPER", STRING, 4);
   MAKEKEY(keywords[KEY_DONOR],    "DONOR",    STRING, 4);
   MAKEKEY(keywords[KEY_ACCEPTOR], "ACCEPTOR", STRING, 1);
   MAKEKEY(keywords[KEY_BUILD],    "BUILD",    STRING, 9);

   /* Parse the file                                                    */
   while(fgets(buffer,MAXBUFF,fp))
   {
      switch(parse(buffer,NKEYS,keywords,RealParam,StrParam))
      {
      case PARSE_ERRC:
         sprintf(gError,"Unrecognised keyword `%s'\n",buffer);
         StoreError("ReadRTop()",gError);
         break;
      case PARSE_ERRP:
         sprintf(gError,"Error in parameters `%s'\n",buffer);
         StoreError("ReadRTop()",gError);
         break;
      case KEY_RESIDUE:
         sscanf(StrParam[1],"%lf",&temp);
         if(!StoreResidueRTop(StrParam[0], temp))
         {
            sprintf(gError,"Too many residue types in %s",filename);
            StoreError("ReadRTop()",gError);
            return(FALSE);
         }
         break;
      case KEY_ATOM:
         sscanf(StrParam[2],"%lf",&temp);
         if(!StoreAtomRTop(StrParam[0], StrParam[1], 
                           temp, StrParam[3]))
         {
            StoreError("ReadRTop()","No memory for atom topology");
            return(FALSE);
         }
         break;
      case KEY_BOND:
         if(!StoreBondRTop(StrParam[0], StrParam[1]))
         {
            StoreError("ReadRTop()","No memory for bond topology");
            return(FALSE);
         }
         break;
      case KEY_ANGLE:
         if(!StoreAngleRTop(StrParam[0], StrParam[1], StrParam[2]))
         {
            StoreError("ReadRTop()","No memory for angle topology");
            return(FALSE);
         }
         break;
      case KEY_TORSION:
         if(!StoreTorsionRTop(StrParam[0], StrParam[1], 
                              StrParam[2], StrParam[3]))
         {
            StoreError("ReadRTop()","No memory for torsion topology");
            return(FALSE);
         }
         break;
      case KEY_IMPROPER:
         if(!StoreImproperRTop(StrParam[0], StrParam[1], 
                               StrParam[2], StrParam[3]))
         {
            StoreError("ReadRTop()","No memory for improper topology");
            return(FALSE);
         }
         break;
      case KEY_DONOR:
         if(!StoreDonorRTop(StrParam[0], StrParam[1], 
                            StrParam[2], StrParam[3]))
         {
            StoreError("ReadRTop()","No memory for donor topology");
            return(FALSE);
         }
         break;
      case KEY_ACCEPTOR:
         if(!StoreAcceptorRTop(StrParam[0]))
         {
            StoreError("ReadRTop()","No memory for acceptor topology");
            return(FALSE);
         }
         break;
      case KEY_BUILD:
         break;
      }
   }

   /* The counter has been used as a pointer into the array rather than
      the actual number of RESIDUE records, so we need to increment.
   */
   gNumResRTop++;
   
   /* Free memory for return strings                                    */
   for(i=0; i<MAXSTRPARAM; i++)
      free(StrParam[i]);

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreResidueRTop(char *res, REAL charge)
   ----------------------------------------------------
   Increments the residue topology count, sets all topology linked list
   pointers to NULL and copies in the residue name and charge.

   30.08.94 Original    By: ACRM
*/
static BOOL StoreResidueRTop(char *res, REAL charge)
{
   if(++gNumResRTop > MAXRESTYPE)
      return(FALSE);

   /* Clear the residue topology entry                                  */
   gRTop[gNumResRTop].AtomTop     = NULL;
   gRTop[gNumResRTop].BondTop     = NULL;
   gRTop[gNumResRTop].AngleTop    = NULL;
   gRTop[gNumResRTop].TorsionTop  = NULL;
   gRTop[gNumResRTop].ImproperTop = NULL;
   gRTop[gNumResRTop].DonorTop    = NULL;
   gRTop[gNumResRTop].AcceptorTop = NULL;

   strcpy(gRTop[gNumResRTop].name, res);
   padterm(gRTop[gNumResRTop].name, 4);

   gRTop[gNumResRTop].charge = charge;

   return(TRUE);
}


/************************************************************************/
/*>static BOOL StoreAtomRTop(char *atom, char *type, REAL charge,
                             char *exclusions)
   --------------------------------------------------------------
   Creates or adds to linked list of atom topology

   30.08.94 Original    By: ACRM
   06.02.03 Fixed for new version of GetWord()
*/
static BOOL StoreAtomRTop(char *atom, char *type, REAL charge,
                          char *exclusions)
{
   ATOMTOP *p;
   char    *excl,
           exclatom[8];
   
   /* Allocate space in linked list for this atom type                  */
   if(gRTop[gNumResRTop].AtomTop == NULL)
   {
      INIT(gRTop[gNumResRTop].AtomTop, ATOMTOP);
      p = gRTop[gNumResRTop].AtomTop;
   }
   else
   {
      p = gRTop[gNumResRTop].AtomTop;
      LAST(p);
      ALLOCNEXT(p, ATOMTOP);
   }

   if(p==NULL)
      return(FALSE);
   
   strcpy(p->atom, atom);
   padterm(p->atom,4);
   strcpy(p->type, type);
   padterm(p->type,4);
   p->charge = charge;

   excl = exclusions;
   p->nexcl=0;
   do 
   {
      excl = GetWord(excl, exclatom, 8);
      padterm(exclatom, 4);
      strcpy(p->excl[p->nexcl], exclatom);

      (p->nexcl)++;
   }  while((excl != NULL) && (p->nexcl<MAXEXCL));

   return(TRUE);
}



/************************************************************************/
/*>static BOOL StoreBondRTop(char *atom1, char *atom2)
   ---------------------------------------------------
   Creates or adds to linked list of bond topology

   30.08.94 Original    By: ACRM
*/
static BOOL StoreBondRTop(char *atom1, char *atom2)
{
   BONDTOP *p;
   
   /* Allocate space in linked list for this bond type                  */
   if(gRTop[gNumResRTop].BondTop == NULL)
   {
      INIT(gRTop[gNumResRTop].BondTop, BONDTOP);
      p = gRTop[gNumResRTop].BondTop;
   }
   else
   {
      p = gRTop[gNumResRTop].BondTop;
      LAST(p);
      ALLOCNEXT(p, BONDTOP);
   }

   if(p==NULL)
      return(FALSE);
   
   strcpy(p->atom1, atom1);
   padterm(p->atom1,4);
   strcpy(p->atom2, atom2);
   padterm(p->atom2,4);

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreAngleRTop(char *atom1, char *atom2, char *atom3)
   -----------------------------------------------------------------
   Creates or adds to linked list of angle topology

   30.08.94 Original    By: ACRM
*/
static BOOL StoreAngleRTop(char *atom1, char *atom2, char *atom3)
{
   ANGLETOP *p;
   
   /* Allocate space in linked list for this angle type                 */
   if(gRTop[gNumResRTop].AngleTop == NULL)
   {
      INIT(gRTop[gNumResRTop].AngleTop, ANGLETOP);
      p = gRTop[gNumResRTop].AngleTop;
   }
   else
   {
      p = gRTop[gNumResRTop].AngleTop;
      LAST(p);
      ALLOCNEXT(p, ANGLETOP);
   }

   if(p==NULL)
      return(FALSE);
   
   strcpy(p->atom1, atom1);
   padterm(p->atom1,4);
   strcpy(p->atom2, atom2);
   padterm(p->atom2,4);
   strcpy(p->atom3, atom3);
   padterm(p->atom3,4);

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreTorsionRTop(char *atom1, char *atom2, 
                                char *atom3, char *atom4)
   ------------------------------------------------------
   Creates or adds to linked list of torsion topology

   30.08.94 Original    By: ACRM
*/
static BOOL StoreTorsionRTop(char *atom1, char *atom2, 
                             char *atom3, char *atom4)
{
   TORSIONTOP *p;
   
   /* Allocate space in linked list for this torsion type               */
   if(gRTop[gNumResRTop].TorsionTop == NULL)
   {
      INIT(gRTop[gNumResRTop].TorsionTop, TORSIONTOP);
      p = gRTop[gNumResRTop].TorsionTop;
   }
   else
   {
      p = gRTop[gNumResRTop].TorsionTop;
      LAST(p);
      ALLOCNEXT(p, TORSIONTOP);
   }

   if(p==NULL)
      return(FALSE);
   
   strcpy(p->atom1, atom1);
   padterm(p->atom1,4);
   strcpy(p->atom2, atom2);
   padterm(p->atom2,4);
   strcpy(p->atom3, atom3);
   padterm(p->atom3,4);
   strcpy(p->atom4, atom4);
   padterm(p->atom4,4);

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreImproperRTop(char *atom1, char *atom2, 
                                 char *atom3, char *atom4)
   -------------------------------------------------------
   Creates or adds to linked list of improper topology

   30.08.94 Original    By: ACRM
*/
static BOOL StoreImproperRTop(char *atom1, char *atom2, 
                              char *atom3, char *atom4)
{
   IMPROPERTOP *p;
   
   /* Allocate space in linked list for this improper type              */
   if(gRTop[gNumResRTop].ImproperTop == NULL)
   {
      INIT(gRTop[gNumResRTop].ImproperTop, IMPROPERTOP);
      p = gRTop[gNumResRTop].ImproperTop;
   }
   else
   {
      p = gRTop[gNumResRTop].ImproperTop;
      LAST(p);
      ALLOCNEXT(p, IMPROPERTOP);
   }

   if(p==NULL)
      return(FALSE);
   
   strcpy(p->atom1, atom1);
   padterm(p->atom1,4);
   strcpy(p->atom2, atom2);
   padterm(p->atom2,4);
   strcpy(p->atom3, atom3);
   padterm(p->atom3,4);
   strcpy(p->atom4, atom4);
   padterm(p->atom4,4);

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreDonorRTop(char *atom1, char *atom2, 
                              char *atom3, char *atom4)
   ----------------------------------------------------
   Creates or adds to linked list of donor topology

   30.08.94 Original    By: ACRM
*/
static BOOL StoreDonorRTop(char *atom1, char *atom2, 
                           char *atom3, char *atom4)
{
   DONORTOP *p;
   
   /* Allocate space in linked list for this donor type                 */
   if(gRTop[gNumResRTop].DonorTop == NULL)
   {
      INIT(gRTop[gNumResRTop].DonorTop, DONORTOP);
      p = gRTop[gNumResRTop].DonorTop;
   }
   else
   {
      p = gRTop[gNumResRTop].DonorTop;
      LAST(p);
      ALLOCNEXT(p, DONORTOP);
   }

   if(p==NULL)
      return(FALSE);
   
   strcpy(p->atom1, atom1);
   padterm(p->atom1,4);
   strcpy(p->atom2, atom2);
   padterm(p->atom2,4);
   strcpy(p->atom3, atom3);
   padterm(p->atom3,4);
   strcpy(p->atom4, atom4);
   padterm(p->atom4,4);

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreAcceptorRTop(char *atom)
   -----------------------------------------
   Creates or adds to linked list of acceptor topology

   30.08.94 Original    By: ACRM
*/
static BOOL StoreAcceptorRTop(char *atom)
{
   ACCEPTORTOP *p;
   
   /* Allocate space in linked list for this Acceptor type              */
   if(gRTop[gNumResRTop].AcceptorTop == NULL)
   {
      INIT(gRTop[gNumResRTop].AcceptorTop, ACCEPTORTOP);
      p = gRTop[gNumResRTop].AcceptorTop;
   }
   else
   {
      p = gRTop[gNumResRTop].AcceptorTop;
      LAST(p);
      ALLOCNEXT(p, ACCEPTORTOP);
   }

   if(p==NULL)
      return(FALSE);
   
   strcpy(p->atom, atom);
   padterm(p->atom,4);

   return(TRUE);
}

/************************************************************************/
/*>void PrintRTop(FILE *out)
   -------------------------
   Prints the contents of the residue topology arrays.

   30.08.94 Original    By: ACRM
*/
void PrintRTop(FILE *out)
{
   int i,
       j;
   ATOMTOP     *pa;
   BONDTOP     *pb;
   ANGLETOP    *pan;
   TORSIONTOP  *pt;
   IMPROPERTOP *pi;
   ACCEPTORTOP *pac;
   DONORTOP    *pd;
   
   for(i=0; i<gNumResRTop; i++)
   {
      fprintf(out,"\nResidue: %4s Charge: %6.3f\n========================\
====\n", gRTop[i].name, gRTop[i].charge);

      fprintf(out,"Atoms:\n");
      for(pa=gRTop[i].AtomTop; pa!=NULL; NEXT(pa))
      {
         fprintf(out,"   %4s Charge: %6.3f Type: %4s ",
                 pa->atom,pa->charge,pa->type);
         for(j=0; j<pa->nexcl; j++)
            fprintf(out,"%4s ",pa->excl[j]);
         fprintf(out,"\n");
      }

      fprintf(out,"Bonds:\n");
      for(pb=gRTop[i].BondTop; pb!=NULL; NEXT(pb))
         fprintf(out,"   %4s %4s\n",pb->atom1,pb->atom2);
         
      fprintf(out,"Angles:\n");
      for(pan=gRTop[i].AngleTop; pan!=NULL; NEXT(pan))
         fprintf(out,"   %4s %4s %4s\n",pan->atom1,pan->atom2,pan->atom3);
         
      fprintf(out,"Torsions:\n");
      for(pt=gRTop[i].TorsionTop; pt!=NULL; NEXT(pt))
         fprintf(out,"   %4s %4s %4s %4s\n",pt->atom1,pt->atom2,pt->atom3,
                 pt->atom4);
         
      fprintf(out,"Impropers:\n");
      for(pi=gRTop[i].ImproperTop; pi!=NULL; NEXT(pi))
         fprintf(out,"   %4s %4s %4s %4s\n",pi->atom1,pi->atom2,pi->atom3,
                 pi->atom4);
         
      fprintf(out,"Donors:\n");
      for(pd=gRTop[i].DonorTop; pd!=NULL; NEXT(pd))
         fprintf(out,"   %4s %4s %4s %4s\n",pd->atom1,pd->atom2,pd->atom3,
                 pd->atom4);
         
      fprintf(out,"Acceptors:\n");
      for(pac=gRTop[i].AcceptorTop; pac!=NULL; NEXT(pac))
         fprintf(out,"   %4s\n",pac->atom);
   }
}

/************************************************************************/
/*>void FreeRTop(void)
   -------------------
   Frees memory allocated for the residue topology data

   30.08.94 Original    By: ACRM
*/
void FreeRTop(void)
{
   int i;
   
   for(i=0; i<gNumResRTop; i++)
   {
      FREELIST(gRTop[i].AtomTop,ATOMTOP);

      FREELIST(gRTop[i].BondTop,BONDTOP);

      FREELIST(gRTop[i].AngleTop,ANGLETOP);

      FREELIST(gRTop[i].TorsionTop,TORSIONTOP);

      FREELIST(gRTop[i].ImproperTop,IMPROPERTOP);

      FREELIST(gRTop[i].DonorTop,DONORTOP);

      FREELIST(gRTop[i].AcceptorTop,ACCEPTORTOP);
   }

   gNumResRTop = 0;
}

/************************************************************************/
/*>int FindRTop(char *resnam)
   --------------------------
   Find the offset into the gRTop residue topology array for a given
   residue name.

   30.08.94 Original    By: ACRM
*/
int FindRTop(char *resnam)
{
   int i;
   
   for(i=0; i<gNumResRTop; i++)
   {
      if(!strncmp(resnam,gRTop[i].name,4))
         return(i);
   }
   
   return(-1);
}

   
   
   

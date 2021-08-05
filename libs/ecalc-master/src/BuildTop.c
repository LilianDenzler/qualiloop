/*************************************************************************

   Program:    ECalc
   File:       BuildTop.c
   
   Version:    V1.5.1
   Date:       07.01.21
   Function:   Build topology from a sequence string
   
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1   30.08.94 Original
   V1.0   30.09.94 First release version
   V1.1   12.10.94 Changes to energy.c
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Skipped
   V1.5.1 07.01.21 Skipped

*************************************************************************/
/* Includes
*/
#define BUILDTOP_MAIN

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/parse.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "bioplib/seq.h"

#include "rtop.h"
#include "ecalc.h"


/************************************************************************/
/* Defines and macros
*/
#define MAXSSSQ ((REAL)5.0)

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
#include "BuildTop.p"

/************************************************************************/
/*>TOPOLOGY *BuildTop(char *seq)
   -----------------------------
   Given a sequence, builds the topology for the whole sequence.
   The input sequence MUST end with a *

   30.08.94 Original    By: ACRM
   14.09.94 Improved error messages
*/
TOPOLOGY *BuildTop(char *seq)
{
   TOPOLOGY *topol;
   int      SeqLen;

   /* Allocate a TOPOLOGY structure to store the data                   */
   if((topol=(TOPOLOGY *)malloc(sizeof(TOPOLOGY)))==NULL)
   {
      StoreError("BuildTop()","No memory for topology structure");
      return(NULL);
   }

   /* Set all pointers to NULL                                          */
   topol->atoms     = NULL;
   topol->bonds     = NULL;
   topol->angles    = NULL;
   topol->torsions  = NULL;
   topol->impropers = NULL;
   topol->donors    = NULL;
   topol->acceptors = NULL;
   topol->ResStart  = NULL;
   topol->sequence  = NULL;

   /* Set all counts to zero                                            */
   topol->NAtoms     = 0;
   topol->NBonds     = 0;
   topol->NAngles    = 0;
   topol->NTorsions  = 0;
   topol->NImpropers = 0;
   topol->NDonors    = 0;
   topol->NAcceptors = 0;
   topol->NRes       = 0;

   topol->Disulphide = FALSE;

   /* Modify the sequence string to add ^ residue to represent 
      NTER (* represents CTER) and place in topology structure
   */
   SeqLen = strlen(seq);
   if((topol->sequence = (char *)malloc((SeqLen+2)*sizeof(char)))==NULL)
   {
      StoreError("BuildTop()","No memory to store sequence");
      DeleteTopology(topol);
      return(NULL);
   }
   topol->sequence[0] = '^';
   strcpy(topol->sequence+1, seq);

   /* Build the atoms array                                             */
   if(BuildAtoms(topol, topol->sequence))
   {
      /* Patch parameters and atom types for the termini                */
      if(PatchTermini(topol))
      {
         /* Build the linked lists used to store all other topology 
            information
         */
         if(BuildRest(topol, topol->sequence))
         {
            /* Patch names of terminal atoms modified by PatchTermini() */
            if(PatchTerminalNames(topol))
            {
               if(PatchNTerDonor(topol))
                  return(topol);
            }
         }
         else
         {
            StoreError("BuildTop()","Unable to build topology");
         }
      }
      else
      {
         StoreError("BuildTop()","Unable to patch termini");
      }
   }
   else
   {
      StoreError("BuildTop()","Unable to store atom topology data");
   }

   DeleteTopology(topol);
   return(NULL);
}

/************************************************************************/
/*>BOOL BuildAtoms(TOPOLOGY *topol, char *sequence)
   ------------------------------------------------
   Build atom data for the sequence into the atom array.

   30.08.94 Original    By: ACRM
   07.09.94 Changed such that topol->atoms is an array of pointers
   14.09.94 Corrected unknown atom message 
   20.09.95 Initialises Mobile flag to TRUE
*/
BOOL BuildAtoms(TOPOLOGY *topol, char *sequence)
{
   int  i, res,
        atom,
        NAtom,
        ResCount;
   char *ch;
   
   /* See how many atoms there are in the sequence                      */
   if((NAtom = NAtomsInSequence(sequence)) == 0)
   {
      StoreError("BuildAtoms()","No atoms in sequence");
      return(FALSE);
   }
   
   /* Now allocate memory for the atom array                            */
   if((topol->atoms=(ATOM **)malloc(NAtom * sizeof(ATOM *)))==NULL)
   {
      StoreError("BuildAtoms()","No memory for topology atom array");
      return(FALSE);
   }

   /* Allocate memory for each ATOM structure                           */
   for(i=0; i<NAtom; i++)
   {
      if((topol->atoms[i]=(ATOM *)malloc(sizeof(ATOM)))==NULL)
      {
         StoreError("BuildAtoms()","No memory for topology atom array");
         return(FALSE);
      }
   }

   /* Now allocate memory for the residue start array                   */
   if((topol->ResStart=(int *)malloc((NAtom+1) * sizeof(int)))==NULL)
   {
      StoreError("BuildAtoms()",
                 "No memory for topology residue offset array");
      return(FALSE);
   }

   /* Zero the atom count                                               */
   topol->NAtoms = 0;
   ResCount      = 0;

   /* Build each residue in the sequence                                */
   for(ch = sequence; *ch; ch++)
   {
      char *resnam;
      ATOMTOP *pa;

      resnam = OneThrTer(*ch);
      
      if((res = FindRTop(resnam)) == (-1))
      {
         sprintf(gError,"%s not in residue topology file",resnam);
         StoreError("BuildAtoms()",gError);
         return(FALSE);
      }
      topol->ResStart[ResCount] = topol->NAtoms;
      for(pa=gRTop[res].AtomTop; pa!=NULL; NEXT(pa))
      {
         strcpy((topol->atoms[topol->NAtoms])->atom, pa->atom);
         strcpy((topol->atoms[topol->NAtoms])->type, pa->type);
         (topol->atoms[topol->NAtoms])->charge = pa->charge;
         if((atom = FindAtomParam(pa->type)) == (-1))
         {
            sprintf(gError,"Unknown atom %s in %s",pa->type,resnam);
            StoreError("BuildAtoms()",gError);
            return(FALSE);
         }
         
         (topol->atoms[topol->NAtoms])->mass   = gAtomParams[atom].mass;
         (topol->atoms[topol->NAtoms])->pol    = gAtomParams[atom].pol;
         (topol->atoms[topol->NAtoms])->NEff   = gAtomParams[atom].NEff;
         (topol->atoms[topol->NAtoms])->vdwr   = gAtomParams[atom].vdwr;
         (topol->atoms[topol->NAtoms])->x      = (REAL)9999.0;
         (topol->atoms[topol->NAtoms])->y      = (REAL)9999.0;
         (topol->atoms[topol->NAtoms])->z      = (REAL)9999.0;
         (topol->atoms[topol->NAtoms])->key    = gAtomParams[atom].key;
         (topol->atoms[topol->NAtoms])->Mobile = TRUE;
         
         (topol->NAtoms)++;
      }
      ResCount++;
   }
   
   /* This is actually beyond the number of residues                    */
   topol->ResStart[ResCount] = topol->NAtoms;

   topol->NRes = ResCount;

   return(TRUE);
}

/************************************************************************/
/*>int NAtomsInSequence(char *sequence)
   ------------------------------------
   Count the atoms in the sequence based on the residue topology file.
   Includes NTER and CTER residues.

   30.08.94 Original    By: ACRM
*/
int NAtomsInSequence(char *sequence)
{
   int  natoms = 0,
        res;
   char *ch;

   for(ch = sequence; *ch; ch++)
   {
      char    *resnam;
      ATOMTOP *pa;

      resnam = OneThrTer(*ch);
      if((res = FindRTop(resnam)) == (-1))
      {
         sprintf(gError,"%s not in residue topology file",resnam);
         StoreError("NAtomsInSequence()",gError);
         return(0);
      }
      for(pa=gRTop[res].AtomTop; pa!=NULL; NEXT(pa))
         natoms++;
   }
   
   return(natoms);
}

/************************************************************************/
/*>BOOL BuildRest(TOPOLOGY *topol, char *sequence)
   -----------------------------------------------
   Builds all topology information other than the atom array.

   30.08.94 Original    By: ACRM
*/
BOOL BuildRest(TOPOLOGY *topol, char *sequence)
{
   char     *ch,
            *resnam;
   int      ResCount=0,
            start,
            stop,
            res;

   for(ch = sequence; *ch; ch++)
   {
      /* Find the atom range for this residue                           */
      start = topol->ResStart[ResCount];
      stop  = topol->ResStart[ResCount+1];

      /* Find the offset into the residue topology array                */
      resnam = OneThrTer(*ch);
      
      if((res = FindRTop(resnam)) == (-1))
      {
         sprintf(gError,"%s not in residue topology",resnam);
         StoreError("BuildRest()",gError);
         return(FALSE);
      }

      if(!BuildExclusions(topol, resnam, res, start, stop))
         return(FALSE);

      if(!BuildBonds(topol, resnam, res, start, stop))
         return(FALSE);

      if(!BuildAngles(topol, resnam, res, start, stop))
         return(FALSE);

      if(!BuildTorsions(topol, resnam, res, start, stop))
         return(FALSE);
      
      if(!BuildImpropers(topol, resnam, res, start, stop))
         return(FALSE);
      
      if(!BuildDonors(topol, resnam, res, start, stop))
         return(FALSE);
      
      if(!BuildAcceptors(topol, resnam, res, start, stop))
         return(FALSE);

      ResCount++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>ATOM *FindAtomInRange(char *atom, TOPOLOGY *topol, int start, int stop)
   -----------------------------------------------------------------------
   Finds the offset into the topology atom array for a specified atom
   name searching between the two specified positions in the array.
   If the atom name starts with + or -, the search begins in the next
   or previous residue

   30.08.94 Original    By: ACRM
   07.09.94 Returns atom pointer rather than offset
            topol->atoms now an array of pointer
   22.05.95 If we're searching for +N and don't find it, then search
            for +OT2 instead.
*/
ATOM *FindAtomInRange(char *atom, TOPOLOGY *topol, int start, int stop)
{
   char atm[8];
   int  i;
   
   switch(*atom)
   {
   case '-':
      strcpy(atm, atom+1);
      padterm(atm,4);
      
      for(i=start-1; i>=0; i--)
      {
         if(!strncmp(atm, topol->atoms[i]->atom, 4))
            return(topol->atoms[i]);
      }
      break;
   case '+':
      strcpy(atm, atom+1);
      padterm(atm,4);
      
      for(i=stop; i<topol->NAtoms; i++)
      {
         if(!strncmp(atm, topol->atoms[i]->atom, 4))
            return(topol->atoms[i]);
      }

      /* If we were searching for +N and didn't find it, look for OT2
         instead.
      */
      if(!strncmp(atm,"N   ",4))
      {
         strcpy(atm,"OT2 ");
         for(i=stop; i<topol->NAtoms; i++)
         {
            if(!strncmp(atm, topol->atoms[i]->atom, 4))
               return(topol->atoms[i]);
         }
      }
      break;
   default:
      strcpy(atm, atom);
      padterm(atm,4);

      for(i=start; i<stop; i++)
      {
         if(!strncmp(atm, topol->atoms[i]->atom, 4))
            return(topol->atoms[i]);
      }
      break;
   }

   return(NULL);
}      

/************************************************************************/
/*>char *OneThrTer(char one)
   -------------------------
   Input:   char  one     One letter code
   Returns: char  *       Three letter code (padded to 4 chars with a 
                          space)

   Converts 1-letter code to 3-letter code (actually as 4 chars).
   This version handles ^ as NTER and * as CTER.

   N.B. The last entry in the static conversion table *must* be UNK.
   
   07.06.93 Original    By: ACRM
   30.08.94 New version for ^ and *
*/
char *OneThrTer(char one)
{
   int j;
   static char sTab1[]    = {'A','C','D','E','F',
                             'G','H','I','K','L',
                             'M','N','P','Q','R',
                             'S','T','V','W','Y',
                             'E','B','Z','^','*',
                             'X',' '
                            };
   static char *sTab3[]   = {"ALA ","CYS ","ASP ","GLU ","PHE ",
                             "GLY ","HIS ","ILE ","LYS ","LEU ",
                             "MET ","ASN ","PRO ","GLN ","ARG ",
                             "SER ","THR ","VAL ","TRP ","TYR ",
                             "PCA ","ASX ","GLX ","NTER","CTER",
                             "UNK ",NULL
                            };

   for(j=0; sTab1[j]!=' '; j++)
      if(sTab1[j] == one) return(sTab3[j]);

   /* Only get here if the one letter code was not found                */
   return(sTab3[j-1]);
}

/************************************************************************/
/*>BOOL PatchTermini(TOPOLOGY *topol)
   -----------------------------------
   Patch the parameters and atom types for the N-ter and C-ter residues

   31.08.94 Original    By: ACRM
   05.09.94 Patches the keys as well
   07.09.94 FindAtomInRange() returns ATOM pointer
            topol->atoms now an array of pointer
*/
BOOL PatchTermini(TOPOLOGY *topol)
{
   int  start, 
        stop,
        i;
   ATOM *atom;
   REAL NTerHCharge;
   
   /* Find start and stop of the true N-terminal residue                */
   start = topol->ResStart[1];
   stop  = topol->ResStart[2];
   
   /* Find the N-terminal nitrogen
      ----------------------------
   */
   if((atom = FindAtomInRange("N   ", topol, start, stop)) == NULL)
   {
      StoreError("PatchTermini()","Unable to patch N-terminal nitrogen");
      return(FALSE);
   }

   /* Patch the atom type and charge                                    */
   strcpy(atom->type,"NH3 ");
   atom->charge += (REAL)0.020;

   /* Find the correct key for this new NH3 type                        */
   for(i=0; i<gNumAtomParams; i++)
   {
      if(!strncmp(gAtomParams[i].atom,"NH3 ",4))
      {
         atom->key = gAtomParams[i].key;
         break;
      }
   }

   /* Find the N-terminal C-alpha
      ---------------------------
   */
   if((atom = FindAtomInRange("CA  ", topol, start, stop)) == NULL)
   {
      StoreError("PatchTermini()","Unable to patch N-terminal C-alpha");
      return(FALSE);
   }

   /* Patch the atom charge                                             */
   atom->charge += (REAL)0.261;

   /* Find the N-terminal hydrogen (if not a proline)
      ----------------------------
   */
   if(topol->sequence[1] != 'P')
   {
      if((atom = FindAtomInRange("H   ", topol, start, stop)) == NULL)
      {
         StoreError("PatchTermini()",
                    "Unable to patch N-terminal hydrogen");
         return(FALSE);
      }

      /* Patch the atom type and parameters                             */
      strcpy(atom->type,"HC  ");

      /* Find the correct key for this new HC type                      */
      for(i=0; i<gNumAtomParams; i++)
      {
         if(!strncmp(gAtomParams[i].atom,"HC  ",4))
         {
            atom->key = gAtomParams[i].key;
            break;
         }
      }
      
      /* Calculate charge                                               */
      NTerHCharge = (atom->charge +
                     topol->atoms[0]->charge    +
                     topol->atoms[1]->charge)/(REAL)3.0;

      /* Patch parameters                                               */
      atom->charge           = NTerHCharge;
      topol->atoms[0]->charge = NTerHCharge;
      topol->atoms[1]->charge = NTerHCharge;
   }
   

   /* Find start and stop of the c-terminal true residue                */
   start = topol->ResStart[topol->NRes - 2];
   stop  = topol->ResStart[topol->NRes - 1];

   /* Find the C-terminal carbon
      --------------------------
   */
   if((atom = FindAtomInRange("C   ", topol, start, stop)) == NULL)
   {
      StoreError("PatchTermini()","Unable to patch C-terminal carbon");
      return(FALSE);
   }
   /* Patch the charge                                                  */
   atom->charge += (REAL)0.030;

   /* Find the C-terminal oxygen
      --------------------------
   */
   if((atom = FindAtomInRange("O   ", topol, start, stop)) == NULL)
   {
      StoreError("PatchTermini()","Unable to patch C-terminal oxygen");
      return(FALSE);
   }
   /* Patch the atom type and charge                                    */
   strcpy(atom->type,"OC  ");
   atom->charge -= (REAL)0.200;

   /* Find the correct key for this new OC type                         */
   for(i=0; i<gNumAtomParams; i++)
   {
      if(!strncmp(gAtomParams[i].atom,"OC  ",4))
      {
         atom->key = gAtomParams[i].key;
         break;
      }
   }

   return(TRUE);
}

/************************************************************************/
/*>BOOL PatchTerminalNames(TOPOLOGY *topol)
   ----------------------------------------
   Patch the names for the terminal atoms whose types and parameters
   were modified by PatchTermini()

   31.08.94 Original    By: ACRM
   07.09.94 FindAtomInRange() returns ATOM pointer
*/
BOOL PatchTerminalNames(TOPOLOGY *topol)
{
   int  start, 
        stop;
   ATOM *atom;
   
   /* Find start and stop of the true N-terminal residue                */
   start = topol->ResStart[1];
   stop  = topol->ResStart[2];
   
   /* Find the N-terminal nitrogen
      ----------------------------
   */
   if((atom = FindAtomInRange("N   ", topol, start, stop)) == NULL)
   {
      StoreError("PatchTerminalNames()",
                 "Unable to patch N-terminal nitrogen");
      return(FALSE);
   }
   strcpy(atom->atom,"NT  ");

   /* Find the N-terminal hydrogen (if not a proline)
      ----------------------------
   */
   if(topol->sequence[1] != 'P')
   {
      if((atom = FindAtomInRange("H   ", topol, start, stop)) == NULL)
      {
         StoreError("PatchTerminalNames()",
                    "Unable to patch N-terminal hydrogen");
         return(FALSE);
      }

      /* Patch the atom type and parameters                             */
      strcpy(atom->atom,"HT3 ");
   }
   
   /* Find start and stop of the c-terminal true residue                */
   start = topol->ResStart[topol->NRes - 2];
   stop  = topol->ResStart[topol->NRes - 1];

   /* Find the C-terminal oxygen
      --------------------------
   */
   if((atom = FindAtomInRange("O   ", topol, start, stop)) == NULL)
   {
      StoreError("PatchTerminalNames()",
                 "Unable to patch C-terminal oxygen");
      return(FALSE);
   }
   /* Patch the atom type and charge                                    */
   strcpy(atom->atom,"OT1 ");

   return(TRUE);
}


/************************************************************************/
/*>BOOL PatchNTerDonor(TOPOLOGY *topol)
   ------------------------------------
   Patches in the donor information for the NTer HT3 atom

   05.09.94 Original    By: ACRM
   07.09.94 FindAtomInRange() returns ATOM pointer
*/
BOOL PatchNTerDonor(TOPOLOGY *topol)
{
   static DONOR *donors;
   ATOM         *atom;
   int          start,
                stop;
   
   if(topol->donors == NULL)
   {
      INIT(topol->donors, DONOR);
      donors = topol->donors;
   }
   else
   {
      donors = topol->donors;
      LAST(donors);
      if(donors->atom1 != NULL)
         ALLOCNEXT(donors, DONOR);
   }

   if(donors==NULL)
   {
      StoreError("PatchNTerDonor()","No memory for donor topology");
      return(FALSE);
   }

   start = topol->ResStart[1];
   stop  = topol->ResStart[2];

   /* Find and store the offsets into the topology for the four
      atoms in this donor
   */
   if((atom = FindAtomInRange("HT3 ", topol, start, stop)) == NULL)
      return(TRUE);
   donors->atom1 = atom;
   if((atom = FindAtomInRange("NT  ", topol, start, stop)) == NULL)
      return(TRUE);
   donors->atom2 = atom;
   if((atom = FindAtomInRange("CA  ", topol, start, stop)) == NULL)
      return(TRUE);
   donors->atom3 = atom;
   if((atom = FindAtomInRange("C   ", topol, start, stop)) == NULL)
      return(TRUE);
   donors->atom4 = atom;
   
   topol->NDonors++;

   return(TRUE);
}


/************************************************************************/
/*>void DeleteTopology(TOPOLOGY *topol)
   ------------------------------------
   Free linked lists and memory allocated for the topology of a chain

   31.08.94 Original    By: ACRM
   08.09.94 Modified for topol->atoms being an array of pointers
*/
void DeleteTopology(TOPOLOGY *topol)
{
   int i;
   
   if(topol != NULL)
   {
      for(i=0; i<topol->NAtoms; i++)
         free(topol->atoms[i]);
      
      if(topol->atoms     != NULL) free(topol->atoms);
      if(topol->bonds     != NULL) FREELIST(topol->bonds, BOND);
      if(topol->angles    != NULL) FREELIST(topol->angles, ANGLE);
      if(topol->torsions  != NULL) FREELIST(topol->torsions, TORSION);
      if(topol->impropers != NULL) FREELIST(topol->impropers, IMPROPER);
      if(topol->donors    != NULL) FREELIST(topol->donors, DONOR);
      if(topol->acceptors != NULL) FREELIST(topol->acceptors, ACCEPTOR);
      free(topol);
   }
   topol = NULL;
}

   

/************************************************************************/
/*>BOOL BuildBonds(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
   -------------------------------------------------------
   Build all bond topology information.

   30.08.94 Original    By: ACRM
   05.09.94 Changed checking on need for allocation
   07.09.94 Changed to store atoms as pointers
*/
BOOL BuildBonds(TOPOLOGY *topol, char *resnam, int res, 
                int start, int stop)
{
   static BOND *bonds = NULL;
   BONDTOP     *bt;
   ATOM        *at1,
               *at2;
   int         BondNum;
   
   /* For each bond in the residue topology                             */
   for(bt=gRTop[res].BondTop; bt!=NULL; NEXT(bt))
   {
      if(topol->bonds==NULL || bonds->atom1 != NULL)
      {
         /* Allocate space in bond linked list                          */
         if(topol->bonds == NULL)
         {
            INIT(topol->bonds, BOND);
            bonds = topol->bonds;
         }
         else
         {
            ALLOCNEXT(bonds, BOND);
         }

         if(bonds==NULL)
         {
            StoreError("BuildBonds()","No memory for bond topology");
            return(FALSE);
         }

         bonds->atom1 = bonds->atom2 = NULL;
      }
      
      
      /* Find and store the offsets into the topology for the two 
         atoms in this bond
      */
      if((at1 = FindAtomInRange(bt->atom1, topol, start, stop)) == NULL)
         continue;
      if((at2 = FindAtomInRange(bt->atom2, topol, start, stop)) == NULL)
         continue;

      bonds->atom1 = at1;
      bonds->atom2 = at2;
      
      /* Find the offset into the bond paramaters array for a bond
         between these two atom types
      */
      if((BondNum = FindBondParam(bonds->atom1->type,
                                  bonds->atom2->type)) 
         == (-1))
      {
         sprintf(gError,"No parameter for bond %s - %s",
                 bonds->atom1->type,
                 bonds->atom2->type);
         StoreError("BuildBonds()",gError);
         return(FALSE);
      }
      
      /* Set the force and optimum length parameters from this
         offset
      */
      bonds->force  = gBondParams[BondNum].force;
      bonds->OptLen = gBondParams[BondNum].OptLen;
      
      topol->NBonds++;
   }

   return(TRUE);
}

         
/************************************************************************/
/*>BOOL BuildAngles(TOPOLOGY *topol, char *resnam, int res, 
                    int start, int stop)
   -------------------------------------------------------
   Build all angle topology information.

   30.08.94 Original    By: ACRM
   05.09.94 Changed checking on need for allocation
   07.09.94 Changed to store atoms as pointers
*/
BOOL BuildAngles(TOPOLOGY *topol, char *resnam, int res, 
                 int start, int stop)
{
   static ANGLE *angles = NULL;
   ANGLETOP     *at;
   ATOM         *at1, *at2, *at3;
   int          AngleNum;
   
   /* For each angle in the residue topology                            */
   for(at=gRTop[res].AngleTop; at!=NULL; NEXT(at))
   {
      if(topol->angles==NULL || angles->atom1 != NULL)
      {
         /* Allocate space in angle linked list                         */
         if(topol->angles == NULL)
         {
            INIT(topol->angles, ANGLE);
            angles = topol->angles;
         }
         else
         {
            ALLOCNEXT(angles, ANGLE);
         }

         if(angles==NULL)
         {
            StoreError("BuildAngles()","No memory for angle topology");
            return(FALSE);
         }

         angles->atom1 = angles->atom2 = angles->atom3 = NULL;
      }
      
      /* Find and store the offsets into the topology for the three
         atoms in this angle
      */
      if((at1 = FindAtomInRange(at->atom1, topol, start, stop)) == NULL)
         continue;
      if((at2 = FindAtomInRange(at->atom2, topol, start, stop)) == NULL)
         continue;
      if((at3 = FindAtomInRange(at->atom3, topol, start, stop)) == NULL)
         continue;

      angles->atom1 = at1;
      angles->atom2 = at2;
      angles->atom3 = at3;
      
      /* Find the offset into the angle paramaters array for an angle
         between these three atom types
      */
      if((AngleNum = FindAngleParam(angles->atom1->type,
                                    angles->atom2->type,
                                    angles->atom3->type)) 
         == (-1))
      {
         sprintf(gError,"No parameter for angle %s - %s - %s",
                 angles->atom1->type,
                 angles->atom2->type,
                 angles->atom3->type);
         StoreError("BuildAngles()",gError);
         return(FALSE);
      }
      
      /* Set the force and optimum angle parameters from this
         offset
      */
      angles->force  = gAngleParams[AngleNum].force;
      angles->OptAng = gAngleParams[AngleNum].OptAng;
      
      topol->NAngles++;
   }

   return(TRUE);
}

/************************************************************************/
/*>BOOL BuildTorsions(TOPOLOGY *topol, char *resnam, int res, 
                      int start, int stop)
   ----------------------------------------------------------
   Build all torsion topology information.

   30.08.94 Original    By: ACRM
   05.09.94 Changed checking on need for allocation
   07.09.94 Changed to store atoms as pointers
*/
BOOL BuildTorsions(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
{
   static TORSION *torsions = NULL;
   TORSIONTOP     *tt;
   ATOM           *at1, *at2, *at3, *at4;
   int            TorsionNum;
   
   /* For each torsion in the residue topology                          */
   for(tt=gRTop[res].TorsionTop; tt!=NULL; NEXT(tt))
   {
      if(topol->torsions==NULL || torsions->atom1 != NULL)
      {
         /* Allocate space in torsion linked list                       */
         if(topol->torsions == NULL)
         {
            INIT(topol->torsions, TORSION);
            torsions = topol->torsions;
         }
         else
         {
            ALLOCNEXT(torsions, TORSION);
         }

         if(torsions==NULL)
         {
            StoreError("BuildTorsions()",
                       "No memory for torsion topology");
            return(FALSE);
         }

         torsions->atom1    = torsions->atom2 = 
            torsions->atom3 = torsions->atom4 = NULL;
      }
      
      /* Find and store the offsets into the topology for the four
         atoms in this torsion
      */
      if((at1 = FindAtomInRange(tt->atom1, topol, start, stop)) == NULL)
         continue;
      if((at2 = FindAtomInRange(tt->atom2, topol, start, stop)) == NULL)
         continue;
      if((at3 = FindAtomInRange(tt->atom3, topol, start, stop)) == NULL)
         continue;
      if((at4 = FindAtomInRange(tt->atom4, topol, start, stop)) == NULL)
         continue;

      torsions->atom1 = at1;
      torsions->atom2 = at2;
      torsions->atom3 = at3;
      torsions->atom4 = at4;
      
      /* Find the offset into the torsion paramaters array for a torsion
         between these four atom types
      */
      if((TorsionNum = 
          FindTorsionParam(torsions->atom1->type,
                           torsions->atom2->type,
                           torsions->atom3->type,
                           torsions->atom4->type)) 
         == (-1))
      {
         sprintf(gError,"No parameter for torsion %s - %s - %s - %s",
                 torsions->atom1->type,
                 torsions->atom2->type,
                 torsions->atom3->type,
                 torsions->atom4->type);
         StoreError("BuildTorsions()",gError);
         return(FALSE);
      }
      
      /* Set the force, period and optimum torsion parameters from this
         offset
      */
      torsions->force  = gTorsionParams[TorsionNum].force;
      torsions->period = gTorsionParams[TorsionNum].period;
      torsions->OptTor = gTorsionParams[TorsionNum].OptTor;
      
      topol->NTorsions++;
   }

   return(TRUE);
}

         
/************************************************************************/
/*>BOOL BuildImpropers(TOPOLOGY *topol, char *resnam, int res, 
                       int start, int stop)
   ----------------------------------------------------------
   Build all improper torsion topology information.

   30.08.94 Original    By: ACRM
   05.09.94 Changed checking on need for allocation
   07.09.94 Changed to store atoms as pointers
*/
BOOL BuildImpropers(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
{
   static IMPROPER *impropers = NULL;
   IMPROPERTOP     *it;
   ATOM            *at1, *at2, *at3, *at4;
   int             ImproperNum;
   
   /* For each improper in the residue topology                         */
   for(it=gRTop[res].ImproperTop; it!=NULL; NEXT(it))
   {
      if(topol->impropers==NULL || impropers->atom1 != NULL)
      {
         /* Allocate space in improper linked list                      */
         if(topol->impropers == NULL)
         {
            INIT(topol->impropers, IMPROPER);
            impropers = topol->impropers;
         }
         else
         {
            ALLOCNEXT(impropers, IMPROPER);
         }

         if(impropers==NULL)
         {
            StoreError("BuildImpropers()",
                       "No memory for improper topology");
            return(FALSE);
         }

         impropers->atom1    = impropers->atom2 = 
            impropers->atom3 = impropers->atom4 = NULL;
      }
      
      /* Find and store the offsets into the topology for the four
         atoms in this improper
      */
      if((at1 = FindAtomInRange(it->atom1, topol, start, stop)) == NULL)
         continue;
      if((at2 = FindAtomInRange(it->atom2, topol, start, stop)) == NULL)
         continue;
      if((at3 = FindAtomInRange(it->atom3, topol, start, stop)) == NULL)
         continue;
      if((at4 = FindAtomInRange(it->atom4, topol, start, stop)) == NULL)
         continue;

      impropers->atom1 = at1;
      impropers->atom2 = at2;
      impropers->atom3 = at3;
      impropers->atom4 = at4;
      
      /* Find the offset into the improper paramaters array for an 
         improper between these four atom types
      */
      if((ImproperNum = 
          FindImproperParam(impropers->atom1->type,
                            impropers->atom2->type,
                            impropers->atom3->type,
                            impropers->atom4->type)) 
         == (-1))
      {
         sprintf(gError,"No parameter for improper %s - %s - %s - %s",
                 impropers->atom1->type,
                 impropers->atom2->type,
                 impropers->atom3->type,
                 impropers->atom4->type);
         StoreError("BuildImpropers()",gError);
         return(FALSE);
      }
      
      /* Set the force and optimum length parameters from this
         offset
      */
      impropers->force  = gImproperParams[ImproperNum].force;
      impropers->OptTor = gImproperParams[ImproperNum].OptTor;
      
      topol->NImpropers++;
   }

   return(TRUE);
}

         
/************************************************************************/
/*>BOOL BuildDonors(TOPOLOGY *topol, char *resnam, int res, 
                       int start, int stop)
   ----------------------------------------------------------
   Build all h-bond donor topology information.

   30.08.94 Original    By: ACRM
   05.09.94 Changed checking on need for allocation
   07.09.94 Changed to store atoms as pointers
*/
BOOL BuildDonors(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
{
   static DONOR *donors = NULL;
   DONORTOP     *dt;
   ATOM         *at1, *at2, *at3, *at4;
   
   /* For each donor in the residue topology                         */
   for(dt=gRTop[res].DonorTop; dt!=NULL; NEXT(dt))
   {
      if(topol->donors==NULL || donors->atom1 != NULL)
      {
         /* Allocate space in donor linked list                      */
         if(topol->donors == NULL)
         {
            INIT(topol->donors, DONOR);
            donors = topol->donors;
         }
         else
         {
            ALLOCNEXT(donors, DONOR);
         }

         if(donors==NULL)
         {
            StoreError("BuildDonors()","No memory for donor topology");
            return(FALSE);
         }

         donors->atom1    = donors->atom2 = 
            donors->atom3 = donors->atom4 = NULL;
      }
      
      /* Find and store the offsets into the topology for the four
         atoms in this donor
      */
      if((at1 = FindAtomInRange(dt->atom1, topol, start, stop)) == NULL)
         continue;
      if((at2 = FindAtomInRange(dt->atom2, topol, start, stop)) == NULL)
         continue;
      if((at3 = FindAtomInRange(dt->atom3, topol, start, stop)) == NULL)
         continue;
      if((at4 = FindAtomInRange(dt->atom4, topol, start, stop)) == NULL)
         continue;

      donors->atom1 = at1;
      donors->atom2 = at2;
      donors->atom3 = at3;
      donors->atom4 = at4;
      
      topol->NDonors++;
   }

   return(TRUE);
}

         
/************************************************************************/
/*>BOOL BuildAcceptors(TOPOLOGY *topol, char *resnam, int res, 
                       int start, int stop)
   ----------------------------------------------------------
   Build all h-bond acceptor topology information.

   30.08.94 Original    By: ACRM
   07.09.94 Changed to store atoms as pointers
*/
BOOL BuildAcceptors(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
{
   static ACCEPTOR *acceptors = NULL;
   ACCEPTORTOP     *at;
   ATOM            *atom;
   
   /* For each acceptor in the residue topology                         */
   for(at=gRTop[res].AcceptorTop; at!=NULL; NEXT(at))
   {
      if(topol->acceptors==NULL || acceptors->atom != NULL)
      {
         /* Allocate space in acceptor linked list                      */
         if(topol->acceptors == NULL)
         {
            INIT(topol->acceptors, ACCEPTOR);
            acceptors = topol->acceptors;
         }
         else
         {
            ALLOCNEXT(acceptors, ACCEPTOR);
         }

         if(acceptors==NULL)
         {
            StoreError("BuildAcceptors()",
                       "No memory for acceptor topology");
            return(FALSE);
         }

         acceptors->atom = NULL;
      }
      
      /* Find and store the offsets into the topology for the four
         atoms in this acceptor
      */
      if((atom = FindAtomInRange(at->atom, topol, start, stop)) == NULL)
         continue;
      acceptors->atom = atom;
      
      topol->NAcceptors++;
   }

   return(TRUE);
}

         
/************************************************************************/
/*>BOOL BuildExclusions(TOPOLOGY *topol, char *resnam, int res, 
                        int start, int stop)
   -------------------------------------------------------
   Build all exclusion topology information.

   30.08.94 Original    By: ACRM
   02.09.94 Modified such that exclusions are atom pointers rather than
            integer offsets
   07.09.94 Changed to store atoms as pointers
*/
BOOL BuildExclusions(TOPOLOGY *topol, char *resnam, int res, 
                     int start, int stop)
{
   int     i,
           NExcl;
   ATOM    *ThisAtom,
           *ExclAtom;
   ATOMTOP *at;
   
   /* Go through each atom for this residue type                        */
   for(at=gRTop[res].AtomTop; at!=NULL; NEXT(at))
   {
      if((ThisAtom = FindAtomInRange(at->atom, topol, start, stop)) 
         != NULL)
      {
         NExcl = 0;
         
         /* Go through the exclusions for this atom                     */
         for(i=0; i<at->nexcl; i++)
         {
            /* Find the excluded atom in the topology                   */
            if((ExclAtom = FindAtomInRange(at->excl[i], topol, 
                                           start, stop)) 
               != NULL)
            {
               /* Add this atom to the exclusion list                   */
               ThisAtom->excl[NExcl++] = ExclAtom;
            }
         }

         ThisAtom->nexcl = NExcl;
      }
   }

   return(TRUE);
}

         
/************************************************************************/
/*>TOPOLOGY *BuildDisulphide(TOPOLOGY *topol1, int res1, 
                             TOPOLOGY *topol2, int res2,
                             MOLECULE *mol, TOPOLOGY *sstopol)
   -----------------------------------------------------------
   Build topology for a disulphide bond between topol1/res1 and
   topol2/res2. If sstopol is NULL, space will be allocated for
   the disulphide topology in the mol->topol array. Otherwise, the
   topology for this disulphide bond will be added to sstopol.

   In any case, the disulphide topology is returned, or NULL if there
   was a fatal error.

   Also adds the atoms involved in the disulphide to the exclusions lists

   N.B. STRICTLY THIS ROUTINE SHOULD PATCH THE PARAMETERS FOR THE 
   ORIGINAL STRUCTURES' BONDS, ETC. AFTER PATCHING THE ATOM TYPE.
   HOWEVER, ALL THE PARAMETERS ARE UNCHANGED IN THE CURRENT PARAMETER
   FILE.

   07.09.94 Original    By: ACRM
*/
TOPOLOGY *BuildDisulphide(TOPOLOGY *topol1, int res1, 
                          TOPOLOGY *topol2, int res2, 
                          MOLECULE *mol, TOPOLOGY *sstopol)
{
   TOPOLOGY *topol = NULL;
   int      start1, stop1,
            start2, stop2,
            BondNum,
            AngleNum,
            TorsionNum;
   ATOM     *s1,    *s2,
            *ca1,   *ca2,
            *cb1,   *cb2;
   static   BOND    *bonds    = NULL;
   static   ANGLE   *angles   = NULL;
   static   TORSION *torsions = NULL;
   static   int     MaxAtoms,
                    AtomNum;

   if(sstopol == NULL)
   {
      /* Set topol to point to next available topology structure        */
      if(mol->NChains >= MAXCHAIN)
      {
         StoreError("BuildDisulphide()",
                    "Too many chains to add disulphide topology");
         return(NULL);
      }

      if((mol->topol[mol->NChains] = 
          (TOPOLOGY *)malloc(sizeof(TOPOLOGY)))==NULL)
      {
         StoreError("BuildDisulphide()",
                    "No memory for disulphide topology");
         return(NULL);
      }
      
      topol = mol->topol[mol->NChains];
      topol->bonds      = NULL;
      topol->atoms      = NULL;
      topol->angles     = NULL;
      topol->torsions   = NULL;
      topol->impropers  = NULL;
      topol->donors     = NULL;
      topol->acceptors  = NULL;
      topol->NRes       = 0;
      topol->NAtoms     = 0;
      topol->NBonds     = 0;
      topol->NAngles    = 0;
      topol->NTorsions  = 0;
      topol->NImpropers = 0;
      topol->NDonors    = 0;
      topol->NAcceptors = 0;

      topol->Disulphide = TRUE;
      
      MaxAtoms          = 0;
      AtomNum           = 0;
      
      (mol->NChains)++;
   }
   else
   {
      topol = sstopol;
   }

   /* Find the atom range for the two residues                          */
   start1 = topol1->ResStart[res1];
   stop1  = topol1->ResStart[res1+1];
   start2 = topol2->ResStart[res2];
   stop2  = topol2->ResStart[res2+1];

   /* Find the involved atoms                                           */
   s1     = FindAtomInRange("SG  ", topol1, start1, stop1);
   cb1    = FindAtomInRange("CB  ", topol1, start1, stop1);
   ca1    = FindAtomInRange("CA  ", topol1, start1, stop1);
   s2     = FindAtomInRange("SG  ", topol2, start2, stop2);
   cb2    = FindAtomInRange("CB  ", topol2, start2, stop2);
   ca2    = FindAtomInRange("CA  ", topol2, start2, stop2);

   if((s1 == NULL) || (cb1 == NULL) || (s2 == NULL) || (cb2 == NULL))
   {
      StoreError("BuildDisulphide()",
                 "Atoms missing from topology data for disulphide bond");
      return(topol);
   }

   /* Patch the exclusions lists for the two sulphurs                   */
   if(!PatchSSExclusions(topol1, s1, cb1, topol2, s2, cb2))
   {
      StoreError("BuildDisulphides()",
                 "Unable to patch exclusions for disulphide bond");
      return(NULL);
   }
   
   /* Add these atoms into the atom list for this topology              */
   if(topol->atoms == NULL)
   {
      if((topol->atoms = (ATOM **)malloc(6 * sizeof(ATOM *)))==NULL)
      {
         StoreError("BuildDisulphide()",
                    "No memory for SS atom topology");
      }
      MaxAtoms = 6;
   }
   else
   {
      MaxAtoms += 6;
      if((topol->atoms = (ATOM **)
          realloc(topol->atoms, MaxAtoms * sizeof(ATOM *)))==NULL)
      {
         StoreError("BuildDisulphide()",
                    "No memory for SS atom topology");
      }
   }

   /* Patch the types for the sulphurs                                  */
   strcpy(s1->type,"S   ");
   strcpy(s2->type,"S   ");
   

   /* Set the atom pointer to point into the main array                 */
   topol->atoms[AtomNum] = s1;
   topol->atoms[AtomNum] = s2;
   topol->atoms[AtomNum] = cb1;
   topol->atoms[AtomNum] = cb2;
   topol->atoms[AtomNum] = ca1;
   topol->atoms[AtomNum] = ca2;


   /* Do the bond
      ===========
   */
   /* Allocate space in bond linked list                                */
   if(topol->bonds == NULL)
   {
      INIT(topol->bonds, BOND);
      bonds = topol->bonds;
   }
   else
   {
      ALLOCNEXT(bonds, BOND);
   }
   
   if(bonds==NULL)
   {
      StoreError("BuildDisulphide()","No memory for SS bond topology");
      return(NULL);
   }
      
   bonds->atom1 = s1;
   bonds->atom2 = s2;
      
   /* Find the offset into the bond paramaters array for a bond
      between these two atom types
   */
   if((BondNum = FindBondParam(s1->type, s2->type)) == (-1))
   {
      sprintf(gError,"No parameter for bond %s - %s",
              s1->type, s2->type);
      StoreError("BuildDisulphide()",gError);
      return(NULL);
   }
      
   /* Set the force and optimum length parameters from this
      offset
   */
   bonds->force  = gBondParams[BondNum].force;
   bonds->OptLen = gBondParams[BondNum].OptLen;
   
   topol->NBonds++;
   
   
   /* Do the first extra angle
      ========================
   */
   /* Allocate space in angle linked list                               */
   if(topol->angles == NULL)
   {
      INIT(topol->angles, ANGLE);
      angles = topol->angles;
   }
   else
   {
      ALLOCNEXT(angles, ANGLE);
   }
   
   if(angles==NULL)
   {
      StoreError("BuildDisulphide()","No memory for SS angle topology");
      return(NULL);
   }
      
   angles->atom1 = cb1;
   angles->atom2 = s1;
   angles->atom3 = s2;
      
   /* Find the offset into the angle paramaters array for a angle
      between these two atom types
   */
   if((AngleNum = FindAngleParam(cb1->type, s1->type, s2->type)) == (-1))
   {
      sprintf(gError,"No parameter for angle %s - %s - %s",
              cb1->type, s1->type, s2->type);
      StoreError("BuildDisulphide()",gError);
      return(NULL);
   }
      
   /* Set the force and optimum length parameters from this
      offset
   */
   angles->force  = gAngleParams[AngleNum].force;
   angles->OptAng = gAngleParams[AngleNum].OptAng;
   
   topol->NAngles++;
   
   /* Do the second extra angle
      =========================
   */
   /* Allocate space in angle linked list                               */
   if(topol->angles == NULL)
   {
      INIT(topol->angles, ANGLE);
      angles = topol->angles;
   }
   else
   {
      ALLOCNEXT(angles, ANGLE);
   }
   
   if(angles==NULL)
   {
      StoreError("BuildDisulphide()","No memory for SS angle topology");
      return(NULL);
   }
      
   angles->atom1 = s1;
   angles->atom2 = s2;
   angles->atom3 = cb2;
      
   /* Find the offset into the angle paramaters array for a angle
      between these two atom types
   */
   if((AngleNum = FindAngleParam(s1->type,s2->type,cb2->type)) == (-1))
   {
      sprintf(gError,"No parameter for angle %s - %s - %s",
              s1->type, s2->type, cb2->type);
      StoreError("BuildDisulphide()",gError);
      return(NULL);
   }
      
   /* Set the force and optimum length parameters from this
      offset
   */
   angles->force  = gAngleParams[AngleNum].force;
   angles->OptAng = gAngleParams[AngleNum].OptAng;
   
   topol->NAngles++;
   
   /* Do the first extra torsion
      ==========================
   */
   /* Allocate space in torsion linked list                             */
   if(topol->torsions == NULL)
   {
      INIT(topol->torsions, TORSION);
      torsions = topol->torsions;
   }
   else
   {
      ALLOCNEXT(torsions, TORSION);
   }
   
   if(torsions==NULL)
   {
      StoreError("BuildDisulphide()","No memory for SS torsion topology");
      return(NULL);
   }
      
   torsions->atom1 = cb1;
   torsions->atom2 = s1;
   torsions->atom3 = s2;
   torsions->atom4 = cb2;
      
   /* Find the offset into the torsion paramaters array for a torsion
      between these two atom types
   */
   if((TorsionNum = FindTorsionParam(cb1->type, s1->type, 
                                     s2->type, cb2->type))
      == (-1))
   {
      sprintf(gError,"No parameter for torsion %s - %s - %s - %s",
              cb1->type,
              s1->type,
              s2->type,
              cb2->type);
      StoreError("BuildDisulphide()",gError);
      return(NULL);
   }
      
   /* Set the force, period and optimum torsion parameters from this
      offset
   */
   torsions->force  = gTorsionParams[TorsionNum].force;
   torsions->period = gTorsionParams[TorsionNum].period;
   torsions->OptTor = gTorsionParams[TorsionNum].OptTor;
   
   topol->NTorsions++;
   
   /* Do the second extra torsion
      ===========================
   */
   /* Allocate space in torsion linked list                             */
   if(topol->torsions == NULL)
   {
      INIT(topol->torsions, TORSION);
      torsions = topol->torsions;
   }
   else
   {
      ALLOCNEXT(torsions, TORSION);
   }
   
   if(torsions==NULL)
   {
      StoreError("BuildDisulphide()","No memory for SS torsion topology");
      return(NULL);
   }
      
   torsions->atom1 = ca1;
   torsions->atom2 = cb1;
   torsions->atom3 = s1;
   torsions->atom4 = s2;
      
   /* Find the offset into the torsion paramaters array for a torsion
      between these two atom types
   */
   if((TorsionNum = FindTorsionParam(ca1->type,
                                     cb1->type,
                                     s1->type,
                                     s2->type))
      == (-1))
   {
      sprintf(gError,"No parameter for torsion %s - %s - %s - %s",
              ca1->type,
              cb1->type,
              s1->type,
              s2->type);
      StoreError("BuildDisulphide()",gError);
      return(NULL);
   }
      
   /* Set the force, period and optimum torsion parameters from this
      offset
   */
   torsions->force  = gTorsionParams[TorsionNum].force;
   torsions->period = gTorsionParams[TorsionNum].period;
   torsions->OptTor = gTorsionParams[TorsionNum].OptTor;
   
   topol->NTorsions++;
   
   /* Do the third extra torsion
      ==========================
   */
   /* Allocate space in torsion linked list                             */
   if(topol->torsions == NULL)
   {
      INIT(topol->torsions, TORSION);
      torsions = topol->torsions;
   }
   else
   {
      ALLOCNEXT(torsions, TORSION);
   }
   
   if(torsions==NULL)
   {
      StoreError("BuildDisulphide()","No memory for SS torsion topology");
      return(NULL);
   }
      
   torsions->atom1 = s1;
   torsions->atom2 = s2;
   torsions->atom3 = cb2;
   torsions->atom4 = ca2;
      
   /* Find the offset into the torsion paramaters array for a torsion
      between these two atom types
   */
   if((TorsionNum = FindTorsionParam(s1->type,
                                     s2->type,
                                     cb2->type,
                                     ca2->type))
      == (-1))
   {
      sprintf(gError,"No parameter for torsion %s - %s - %s - %s",
              s1->type,
              s2->type,
              cb2->type,
              ca2->type);
      StoreError("BuildDisulphide()",gError);
      return(NULL);
   }
      
   /* Set the force, period and optimum torsion parameters from this
      offset
   */
   torsions->force  = gTorsionParams[TorsionNum].force;
   torsions->period = gTorsionParams[TorsionNum].period;
   torsions->OptTor = gTorsionParams[TorsionNum].OptTor;
   
   topol->NTorsions++;
   
   return(topol);
}


/************************************************************************/
/*>int ScanDisulphides(MOLECULE *mol)
   -----------------------------------
   Scan through the complete structure to search for disulphides using
   a simple distance criterion. When a disulphide is found, calls the
   BuildDisulphide() routine to create its topology and to patch the
   parent structure topology.

   Returns number of disulphide found or -1 for error.

   07.09.94 Original    By: ACRM
   08.09.94 Modified for topol->atoms being an array of pointers
   09.09.94 Fixed bug in scanning other chains for Cys-SG
*/
int ScanDisulphides(MOLECULE *mol)
{
   int      i,    j,
            res1, res2,
            AtomCount1,
            AtomCount2,
            NumSS = 0;
   TOPOLOGY *topol1,
            *topol2,
            *sstopol = NULL;
   ATOM     *atom1,
            *atom2;
   

   /* First dispose of any current disulphide topology                  */
   for(i = 0; i < mol->NChains; i++)
   {
      topol1 = mol->topol[i];
      
      if(topol1->Disulphide)
      {
         /* Note: we do *not* free the individual atom entries as these
            are simply references to the main topology arrays
         */
         free(topol1->atoms);
         
         FREELIST(topol1->angles, ANGLE);
         FREELIST(topol1->torsions, TORSION);
         free(topol1);
         mol->NChains = i-1;
         break;
      }
   }
   
   /* For each chain in the structure                                   */
   for(i=0; i<mol->NChains; i++)
   {
      topol1 = mol->topol[i];
      
      /* For each residue in the chain                                  */
      for(res1 = 0; res1<topol1->NRes; res1++)
      {
         /* If it's a Cysteine                                          */
         if(topol1->sequence[res1] == 'C')
         {
            /* For each atom in this residue                            */
            for(AtomCount1 = topol1->ResStart[res1];
                AtomCount1 < topol1->ResStart[res1+1];
                AtomCount1++)
            {
               atom1 = topol1->atoms[AtomCount1];

               /* If got a CYS SG                                       */
               if(!strncmp(atom1->atom, "SG  ",4))
               {
                  /* For each other residue in the chain                */
                  for(res2 = res1+1; res2<topol1->NRes; res2++)
                  {
                     /* If it's a Cysteine                              */
                     if(topol1->sequence[res2] == 'C')
                     {
                        /* For each atom in this residue                */
                        for(AtomCount2 = topol1->ResStart[res2];
                            AtomCount2 < topol1->ResStart[res2+1];
                            AtomCount2++)
                        {
                           atom2 = topol1->atoms[AtomCount2];

                           /* If got a CYS SG                           */
                           if(!strncmp(atom2->atom, "SG  ",4))
                           {
                              /* See if it's in range                   */
                              if(DISTSQ(atom1, atom2) < MAXSSSQ)
                              {
                                 sstopol = 
                                    BuildDisulphide(topol1, res1, 
                                                    topol1, res2,
                                                    mol, sstopol);
                                 if(sstopol == NULL)
                                    return(-1);

                                 NumSS++;
                              }
                           }
                        }
                     }
                  }
                  
                  /* For each remaining chain                           */
                  for(j=i+1; j<mol->NChains; j++)
                  {
                     topol2 = mol->topol[j];
      
                     /* For each residue in the chain                   */
                     for(res2 = 0; res2<topol2->NRes; res2++)
                     {
                        /* If it's a Cysteine                           */
                        if(topol2->sequence[res2] == 'C')
                        {
                           /* For each atom in this residue             */
                           for(AtomCount2 = topol2->ResStart[res2];
                               AtomCount2 < topol2->ResStart[res2+1];
                               AtomCount2++)
                           {
                              atom2 = topol2->atoms[AtomCount2];

                              /* If got a CYS SG                        */
                              if(!strncmp(atom2->atom, "SG  ",4))
                              {
                                 /* See if it's in range                */
                                 if(DISTSQ(atom1, atom2) < MAXSSSQ)
                                 {
                                    sstopol = 
                                       BuildDisulphide(topol1, res1, 
                                                       topol2, res2,
                                                       mol, sstopol);
                                    if(sstopol == NULL)
                                       return(-1);
                                    NumSS++;
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return(NumSS);
}

/************************************************************************/
/*>BOOL PatchSSExclusions(TOPOLOGY *topol1, ATOM *AtomS1, ATOM *AtomCB1,
                          TOPOLOGY *topol2, ATOM *AtomS2, ATOM *AtomCB2)
   ---------------------------------------------------------------------
   Patch the exclusions associated with a disulphide bond

   09.09.94 Original    By: ACRM
*/
BOOL PatchSSExclusions(TOPOLOGY *topol1, ATOM *AtomS1, ATOM *AtomCB1,
                       TOPOLOGY *topol2, ATOM *AtomS2, ATOM *AtomCB2)
{
   if((AtomS1->nexcl >= MAXEXCL) || (AtomCB1->nexcl >= MAXEXCL) ||
      (AtomS2->nexcl >= MAXEXCL) || (AtomCB2->nexcl >= MAXEXCL))
   {
      StoreError("PatchSSExclusions() [Increase MAXEXCL]",
                 "Too many exclusions");
      return(FALSE);
   }
   
   /* Exclude S2 from S1                                                */
   AtomS1->excl[AtomS1->nexcl] = AtomS2;
   if(++(AtomS1->nexcl) >= MAXEXCL)
   {
      StoreError("PatchSSExclusions() [Increase MAXEXCL]",
                 "Too many exclusions");
      return(FALSE);
   }

   /* Exclude S2 from CB1                                               */
   AtomCB1->excl[AtomCB1->nexcl] = AtomS2;
   if(++(AtomCB1->nexcl) >= MAXEXCL)
   {
      StoreError("PatchSSExclusions() [Increase MAXEXCL]",
                 "Too many exclusions");
      return(FALSE);
   }

   /* Exclude S1 from CB2                                               */
   AtomCB2->excl[AtomCB2->nexcl] = AtomS1;
   if(++(AtomCB2->nexcl) >= MAXEXCL)
   {
      StoreError("PatchSSExclusions() [Increase MAXEXCL]",
                 "Too many exclusions");
      return(FALSE);
   }

   return(TRUE);
}


/************************************************************************/
/*>int GetResFromAtomNum(TOPOLOGY *topol, int atomnum)
   ---------------------------------------------------
   Gets the residue number (counting from 0) for an atom number in a
   given topology. Returns -1 if atom num out of range

   13.09.94 Original    By: ACRM
*/
int GetResFromAtomNum(TOPOLOGY *topol, int atomnum)
{
   int i;
   
   for(i=0; i<topol->NRes; i++)
   {
      if(atomnum < topol->ResStart[i])
      {
         return(i-1);
      }
   }
   
   return(-1);
}

/************************************************************************/
/*>BOOL UserDisulphide(MOLECULE *mol, char *resspec1, char *resspec2)
   ------------------------------------------------------------------
   Handle a DISULPHIDE command given by the user to patch a disulphide 
   manually

   13.09.94 Original    By: ACRM
*/
BOOL UserDisulphide(MOLECULE *mol, char *resspec1, char *resspec2)
{
   char     chain1,
            chain2,
            insert1,
            insert2;
   int      i,
            resnum1,
            resnum2;
   TOPOLOGY *topol1  = NULL,
            *topol2  = NULL;
   static TOPOLOGY *sstopol = NULL;

   /* Split up the residue specifications                               */
   ParseResSpec(resspec1, &chain1, &resnum1, &insert1);
   ParseResSpec(resspec2, &chain2, &resnum2, &insert2);

   /* Decrement the residue numbers since we will count from 0          */
   resnum1--;
   resnum2--;

   /* Check that consecutive numbering has been used                    */
   if(insert1 != ' ' || insert2 != ' ')
   {
      StoreError("UserDisulphide()",
                 "Consecutive numbering (no insertions) must be used in \
DISULPHIDE commands");
      return(FALSE);
   }

   /* Find the topologies for the two chains                            */
   for(i=0; i<mol->NChains; i++)
   {
      if(!(mol->topol[i]->Disulphide))
      {
         if(topol1==NULL && mol->topol[i]->ChainName == chain1)
            topol1 = mol->topol[i];
         if(topol2==NULL && mol->topol[i]->ChainName == chain2)
            topol2 = mol->topol[i];
      }
   }
   
   /* Check we found both chains                                        */
   if(topol1==NULL)
   {
      sprintf(gError,"Invalid chain (%c) in DISULPHIDE command",chain1);
      StoreError("UserDisulphide()",gError);
      return(FALSE);
   }
   if(topol2==NULL)
   {
      sprintf(gError,"Invalid chain (%c) in DISULPHIDE command",chain2);
      StoreError("UserDisulphide()",gError);
      return(FALSE);
   }

   /* Check that these residues are cysteines                           */
   if(topol1->sequence[resnum1] != 'C')
   {
      sprintf(gError,"First residue in DISULPHIDE command (%s) is not a \
CYS", resspec1);
      StoreError("UserDisulphide()",gError);
      return(FALSE);
   }
   if(topol2->sequence[resnum2] != 'C')
   {
      sprintf(gError,"Second residue in DISULPHIDE command (%s) is not a \
CYS", resspec2);
      StoreError("UserDisulphide()",gError);
      return(FALSE);
   }

   /* Patch in the disulphide                                           */
   if((sstopol = BuildDisulphide(topol1, resnum1, topol2, resnum2, 
                                 mol, sstopol)) == NULL)
   {
      StoreError("UserDisulphide()","Unable to build disulphide");
      return(FALSE);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL PatchZones(MOLECULE *mol, ZONE *zones)
   -------------------------------------------
   Runs through a set of zones specified as residue specs (e.g. L24)
   and adds topology and atom pointers

   21.09.94 Extracted from TrimTopology()  By: ACRM
   22.09.94 Corrected error message
   31.01.95 Initialised topol for gcc
   19.05.95 Corrected call to StoreError() to specify PatchZones() not
            TrimTopology()
*/
BOOL PatchZones(MOLECULE *mol, ZONE *zones)
{
   ZONE     *z;
   char     chain1,  chain2,
            insert1, insert2;
   int      resnum1, resnum2,
            count;
   TOPOLOGY *topol = NULL;

   if(zones==NULL)
      return(TRUE);

   /* Patch the ZONE linked list to get atom number ranges              */
   for(z=zones; z!=NULL; NEXT(z))
   {
      /* Split up the residue specifications                            */
      ParseResSpec(z->res1, &chain1, &resnum1, &insert1);
      ParseResSpec(z->res2, &chain2, &resnum2, &insert2);

      /* Decrement the residue numbers since we will count from 0       */
      resnum1--;
      resnum2--;

      /* Find the topology for the chain (we have already checked the
         chains are identical when storing these data)
      */
      for(count=0; count<mol->NChains; count++)
      {
         if(!(mol->topol[count]->Disulphide))
         {
            if(mol->topol[count]->ChainName == chain1)
            {
               topol = mol->topol[count];
               break;
            }
         }
      }
   
      /* Check we found the chain                                       */
      if(topol==NULL)
      {
         sprintf(gError,"Invalid chain specification (%c)",chain1);
         StoreError("PatchZones()",gError);
         return(FALSE);
      }

      /* Store the topology in the ZONE linked list                     */
      z->topol = topol;

      /* Find the first and following atom numbers in the range         */
      z->atom1 = topol->ResStart[resnum1];
      z->atom2 = topol->ResStart[resnum2+1];
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL TrimTopology(MOLECULE *mol, ZONE *zones)
   ---------------------------------------------
   Removes all bonded topology information from the topology lists
   which include no atoms in the specified zones.

   Also sets the Mobile flag for each atom in a zone.

   16.09.94 Original    By: ACRM
   21.09.94 Removed the zone patching code which is now in PatchZones()
            and is called after the structure has been read.
   29.09.94 Removed unused variables
   19.05.95 Runs through the atom lists and sets Mobile flag for all atoms
            in the zones
*/
BOOL TrimTopology(MOLECULE *mol, ZONE *zones)
{
   int       count,
             AtomNum;
   TOPOLOGY  *topol  = NULL;
   BOND      *b, *lastb;
   ANGLE     *a, *lasta;
   TORSION   *t, *lastt;
   IMPROPER  *i, *lasti;

   if(zones == NULL)
      return(TRUE);

   /* For each topology                                                 */
   for(count=0; count<mol->NChains; count++)
   {
      topol = mol->topol[count];

      /* Set all last pointers to NULL                                  */
      lastb = NULL;
      lasta = NULL;
      lastt = NULL;
      lasti = NULL;

      /* Set all topology counts to 0                                   */
      topol->NBonds     = 0;
      topol->NAngles    = 0;
      topol->NTorsions  = 0;
      topol->NImpropers = 0;

      /* 19.05.95 Set mobile flag for each atom                         */
      for(AtomNum=0; AtomNum<topol->NAtoms; AtomNum++)
      {
         if(topol->atoms[AtomNum])
         {
            if(InZones(topol->atoms[AtomNum], zones))
            {
               (topol->atoms[AtomNum])->Mobile = TRUE;
            }
            else
            {
               (topol->atoms[AtomNum])->Mobile = FALSE;
            }
         }
      }
      
      /* Run through the BOND linked list checking that entries have 
         at least one atom in one of the zones
      */
      for(b=topol->bonds; b!=NULL; )
      {
         /* If neither atom is in a zone, unlink from linked list       */
         if(!InZones(b->atom1, zones) && !InZones(b->atom2, zones))
         {
            if(lastb == NULL)
            {
               topol->bonds = b->next;
               free(b);
               b = topol->bonds;
            }
            else
            {
               lastb->next = b->next;
               free(b);
               b = lastb->next;
            }
         }
         else      /* Keep this one; just update last pointer           */
         {
            lastb = b;
            NEXT(b);
            (topol->NBonds)++;
         }
      }
      
      /* Run through the ANGLE linked list checking that entries have 
         at least one atom in one of the zones
      */
      for(a=topol->angles; a!=NULL; )
      {
         /* If no atom is in a zone, unlink from linked list            */
         if(!InZones(a->atom1, zones) && 
            !InZones(a->atom2, zones) &&
            !InZones(a->atom3, zones))
         {
            if(lasta == NULL)
            {
               topol->angles = a->next;
               free(a);
               a = topol->angles;
            }
            else
            {
               lasta->next = a->next;
               free(a);
               a = lasta->next;
            }
         }
         else      /* Keep this one; just update last pointer           */
         {
            lasta = a;
            NEXT(a);
            (topol->NAngles)++;
         }
      }
         
      /* Run through the TORSION linked list checking that entries have 
         at least one atom in one of the zones
      */
      for(t=topol->torsions; t!=NULL; )
      {
         /* If no atom is in a zone, unlink from linked list            */
         if(!InZones(t->atom1, zones) && 
            !InZones(t->atom2, zones) &&
            !InZones(t->atom3, zones) &&
            !InZones(t->atom4, zones))
         {
            if(lastt == NULL)
            {
               topol->torsions = t->next;
               free(t);
               t = topol->torsions;
            }
            else
            {
               lastt->next = t->next;
               free(t);
               t = lastt->next;
            }
         }
         else      /* Keep this one; just update last pointer           */
         {
            lastt = t;
            NEXT(t);
            (topol->NTorsions)++;
         }
      }
      
      /* Run through the IMPROPER linked list checking that entries have 
         at least one atom in one of the zones
      */
      for(i=topol->impropers; i!=NULL; )
      {
         /* If no atom is in a zone, unlink from linked list            */
         if(!InZones(i->atom1, zones) && 
            !InZones(i->atom2, zones) &&
            !InZones(i->atom3, zones) &&
            !InZones(i->atom4, zones))
         {
            if(lasti == NULL)
            {
               topol->impropers = i->next;
               free(i);
               i = topol->impropers;
            }
            else
            {
               lasti->next = i->next;
               free(i);
               i = lasti->next;
            }
         }
         else      /* Keep this one; just update last pointer           */
         {
            lasti = i;
            NEXT(i);
            (topol->NImpropers)++;
         }
      }
   }  /* Next topology                                                  */
   
   return(TRUE);
}



/************************************************************************/
/*>BOOL InZones(ATOM *atom, ZONE *zones)
   -------------------------------------
   Runs through the ZONE linked list to see if the given ATOM pointer
   is within any of the specified zones. If zones is NULL, allways
   returns TRUE
   The zones structure also contains an `sc' flag. If this is true,
   the return will only be TRUE when atom is in a zone AND is a s/c
   atom (from CG out). This is used when the zones are ignore-zones
   specified with the sidechain qualifier.

   16.09.94 Original    By: ACRM
   21.09.94 Added checking of sc flag in zones
*/
BOOL InZones(ATOM *atom, ZONE *zones)
{
   ZONE *z;
   int  AtomNum;

   if(zones == NULL)
      return(TRUE);

   for(z=zones; z!=NULL; NEXT(z))
   {
      for(AtomNum = z->atom1; AtomNum < z->atom2; AtomNum++)
      {
         if(atom == z->topol->atoms[AtomNum])
         {
            if(!z->sc)       /* Normal behaviour, any atom valid        */
            {
               return(TRUE);
            }
            else
            {
               /* If it's a backbone atom, return FALSE; sidechain TRUE */
               if(!strncmp(atom->atom,"N   ",4) ||
                  !strncmp(atom->atom,"H   ",4) ||
                  !strncmp(atom->atom,"CA  ",4) ||
                  !strncmp(atom->atom,"C   ",4) ||
                  !strncmp(atom->atom,"O   ",4) ||
                  !strncmp(atom->atom,"CB  ",4) ||
                  !strncmp(atom->atom,"NT  ",4) ||
                  !strncmp(atom->atom,"HT1 ",4) ||
                  !strncmp(atom->atom,"HT2 ",4) ||
                  !strncmp(atom->atom,"HT3 ",4) ||
                  !strncmp(atom->atom,"OT1 ",4) ||
                  !strncmp(atom->atom,"OT2 ",4))
               {
                  return(FALSE);
               }
               else
               {
                  return(TRUE);
               }
            }
         }
      }
   }
   
   return(FALSE);
}

      

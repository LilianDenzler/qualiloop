/*************************************************************************

   Program:    ECalc
   File:       ReadStructure.c
   
   Version:    V1.5.2
   Date:       05.03.21
   Function:   Read sequence & coordinates into the topology structure
   
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
   V0.1   01.09.94 Original
   V1.0   30.09.94 First release version
   V1.1   12.10.94 Changes to energy.c
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Allow GLU/PCA to match
                   Added code to remove PRO/HT3
   V1.5.1 07.01.21 Skipped
   V1.5.2 05.03.21 Skipped

*************************************************************************/
/* Includes
*/
#define READSTRUCTURE_MAIN

#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/angle.h"
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"

#include "ecalc.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>MOLECULE *ReadStructure(FILE *fp, BOOL DoSS, BOOL *MissingAtoms, 
                           int *NumSS)
   ----------------------------------------------------------------
   Build a molecule structure given a pointer to a PDB file

   01.09.94 Original    By: ACRM
   02.09.94 Added mol->NAtoms
   05.09.94 Added mol->NDonors
   07.09.94 Added ScanDisulphides()
   13.09.93 Outputs number of disulphides found and only does disulphide
            search if the DoSS flag is set
            Added code to handle a user-supplied disulphide list
   14.09.94 Corrected error message when BuildTop() fails
   16.09.94 Changed name of USERSS to ZONE
            Added call to TrimTopology()
   21.09.94 Added calls to PatchZones(). This was previously done
            within TrimTopology().
   22.09.94 Checks error returns from PatchZones()
   06.02.03 Added checks on allocations before free()s
   07.02.03 Added code to set PRO/HT3 to NULL coords
*/
MOLECULE *ReadStructure(FILE *fp, BOOL DoSS, BOOL *MissingAtoms, 
                        int *NumSS)
{
   PDB      *pdb,
            *p;
   int      natoms,
            i, 
            NSeq;
   char     *sequence,
            *seqs[MAXCHAIN];
   MOLECULE *mol;
   BOOL     error;

   /* Read structure from file into PDB linked list                     */
   if((pdb=ReadPDB(fp, &natoms))==NULL)
   {
      StoreError("ReadStructure()","No atoms read from PDB file");
      return(NULL);
   }

   /* 07.02.03 Set any PRO/HT3 atoms to NULL coords - actually PRO/HT3
      seems to be ignored correctly anyway, but this is a useful
      precaution!
   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->resnam, "PRO ", 4) &&
         !strncmp(p->atnam,  "HT3 ", 4))
      {
         p->x = p->y = p->z = 9999.0;
      }
   }

   /* Allocate the molecule structure                                   */
   if((mol = (MOLECULE *)malloc(sizeof(MOLECULE)))==NULL)
   {
      FREELIST(pdb,PDB);
      StoreError("ReadStructure()","No memory for molecule structure");
      return(NULL);
   }
   mol->NAtoms = 0;
   mol->NDonors = 0;
   
   /* Extract sequence from PDB linked list. N.B. This routine ignores
      NTER and CTER residues
   */
   if((sequence = PDB2Seq(pdb))==NULL)
   {
      FREELIST(pdb, PDB);
      free(mol);
      StoreError("ReadStructure()","No memory for sequence");
      return(NULL);
   }

   /* Split this into separate chains                                   */
   if((NSeq = SplitSeq(sequence, seqs))==0)
   {
      FREELIST(pdb,PDB);
      free(sequence);
      free(mol);
      StoreError("ReadStructure()","No memory for split sequences");
      return(NULL);
   }
   mol->NChains = NSeq;

   /* Build topology for each chain                                     */
   for(i=0; i<NSeq; i++)
   {
      char *StarSeq;
      
      /* Make a copy of the sequence with a * at the end                */
      if((StarSeq=(char *)malloc((strlen(seqs[i])+2)*sizeof(char)))==NULL)
      {
         FREELIST(pdb,PDB);
         free(sequence);
         for(i=0; i<NSeq; i++)
            free(seqs[i]);
         StoreError("ReadStructure()","No memory to copy sequence.");
         return(NULL);
      }
      strcpy(StarSeq,seqs[i]);
      strcat(StarSeq,"*");
      
      /* Actually build the topology                                    */
      /* 06.02.03 Added checks before frees                             */
      if((mol->topol[i] = BuildTop(StarSeq))==NULL)
      {
         if(pdb) FREELIST(pdb,PDB);
         if(sequence) free(sequence);
         for(i=0; i<NSeq; i++)
         {
            if(seqs[i])  free(seqs[i]);
         }
         
         StoreError("ReadStructure()","Unable to build topology.");
         return(NULL);
      }

      /* Free the allocated sequence space                              */
      free(StarSeq);
   }

   /* Store coordinates for each chain                                  */
   error = FALSE;
   for(i=0, p=pdb; i<NSeq && p!=NULL; i++)
      p = StoreCoords(mol->topol[i], p, &error);

   if(error)
   {
      StoreError("ReadStructure()","Unable to store coordinates");
      return(NULL);
   }

   /* Check that all coordinates were read successfully                 */
   for(i=0; i<NSeq; i++)
      if(!CheckCoords(mol->topol[i])) *MissingAtoms = TRUE;

   /* Free the PDB linked list                                          */
   FREELIST(pdb,PDB);

   /* Count the atoms and donors from the mol structure                 */
   for(i=0; i<NSeq; i++)
   {
      if(!(mol->topol[i])->Disulphide)
      {
         mol->NAtoms  += (mol->topol[i])->NAtoms;
         mol->NDonors += (mol->topol[i])->NDonors;
      }
   }

   if(DoSS)     /* Auto disulphide search                               */
   {
      /* Scan the structure for disulphides and add their topology      */
      if((*NumSS = ScanDisulphides(mol)) == (-1))
      {
         StoreError("ReadStructure()",
                    "Failure in scanning for disulphides");
         return(NULL);
      }
   }
   else         /* Any user supplied disulphides                        */
   {
      ZONE *s;
      
      for(s=gUserSS; s!=NULL; NEXT(s))
      {
         if(!UserDisulphide(mol, s->res1, s->res2))
         {
            StoreError("ReadStructure()",
                       "Failed to build user disulphide");
            return(NULL);
         }
      }
   }

   /* Patch the zone and ignore linked lists with the topology and 
      atom pointer information
   */
   if(!PatchZones(mol, gZone))
   {
      StoreError("ReadStructure()",
                 "Error in ZONE command");
      return(NULL);
   }

   if(!PatchZones(mol, gIgnore))
   {
      StoreError("ReadStructure()",
                 "Error in IGNORE command");
      return(NULL);
   }

   /* Now trim the topology down to any user supplied zones             */
   if(!TrimTopology(mol, gZone))
   {
      StoreError("ReadStructure()",
                 "Failed to trim topology to specified zones");
      return(NULL);
   }

   return(mol);
}

/************************************************************************/
/*>PDB *StoreCoords(TOPOLOGY *topol, PDB *pdb, BOOL *error)
   --------------------------------------------------------
   Fills in the coordinates in the topology structure from the pdb
   linked list. Returns a pointer into the linked list after all the
   coordinates which have been used so far

   01.09.94 Original    By: ACRM
   08.09.94 Changed for atoms as array of pointers
   13.09.94 Stores the chain name in the topol structure
   06.02.03 Don't raise an error if it's GLU/PCA
*/
PDB *StoreCoords(TOPOLOGY *topol, PDB *pdb, BOOL *error)
{
   int  i, j,
        AtStart,
        AtStop;
   PDB  *p,
        *pdbstart,
        *pdbstop;
   char chain,
        *resnam;

   pdbstart = pdb;
   chain    = pdb->chain[0];

   topol->ChainName = chain;
   
   /* For each residue in the topology                                  */
   for(i=0; i<topol->NRes; i++)
   {
      /* Check the residue type matches                                 */
      resnam = OneThrTer(topol->sequence[i]);
      if(strncmp(pdbstart->resnam, resnam, 4) &&
         (strncmp(pdbstart->resnam, "PCA ", 4) ||
          strncmp(resnam, "GLU ", 4)))
      {
         sprintf(gError,"Residue mismatch (Topology %s, coords %s)",
                 resnam, pdbstart->resnam);
         StoreError("StoreCoords()",gError);
         *error=TRUE;
         return(NULL);
      }
      
      /* Find the atom range for this residue                           */
      AtStart = topol->ResStart[i];
      AtStop  = topol->ResStart[i+1];

      /* Find the end of this residue in the PDB linked list            */
      pdbstop = FindEndPDB(pdbstart);
      
      /* For each atom in the topology                                  */
      for(j=AtStart; j<AtStop; j++)
      {
         /* Run through the PDB linked list                             */
         for(p=pdbstart; p!=pdbstop; NEXT(p))
         {
            if(p->chain[0] != chain)
            {
               StoreError("StoreCoords()","Chain in PDB has ended too \
early");
               return(NULL);
            }

            if(!strncmp((topol->atoms[j])->atom, p->atnam, 4))
            {
               (topol->atoms[j])->x = p->x;
               (topol->atoms[j])->y = p->y;
               (topol->atoms[j])->z = p->z;
            }
            else if(!strncmp(resnam,"ILE ",4) &&
                    ((!strncmp((topol->atoms[j])->atom,"CD  ",4) &&
                      !strncmp(p->atnam,"CD1 ",4)) ||
                     (!strncmp((topol->atoms[j])->atom,"CD1 ",4) &&
                      !strncmp(p->atnam,"CD  ",4))))
            {
               /* This is an ILE with CD matching CD1                   */
               (topol->atoms[j])->x = p->x;
               (topol->atoms[j])->y = p->y;
               (topol->atoms[j])->z = p->z;
            }
         }
      }

      pdbstart = pdbstop;
   }

   return(pdbstart);
}

/************************************************************************/
/*>BOOL CheckCoords(TOPOLOGY *topol)
   ---------------------------------
   Run through the topology atom records and check they all have
   coordinates set.

   01.09.94 Original   By: ACRM
   08.09.94 Changed for atoms as array of pointers
*/
BOOL CheckCoords(TOPOLOGY *topol)
{
   int i;

   for(i=0; i<topol->NAtoms; i++)
   {
      if((topol->atoms[i]->x > 9998.0) &&
         (topol->atoms[i]->y > 9998.0) &&
         (topol->atoms[i]->z > 9998.0))
         return(FALSE);
   }
   
   return(TRUE);
}


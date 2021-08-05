#include "ProtoTypes.h"
#include "CongenProto.h"
 
/***************************************************************************
 
   Program:
   File:       hadd.c
 
   Version:    V2.3
   Date:
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
 
   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission from
   the author, although it may be given away free with commercial products,
   providing it is made clear that this program is free and that the source
   code is provided with the program.
 
   Description:
   ============
   Routine to add hydrogens to a protein linked list of type PDB.
   The routine allocates space for the new atoms and inserts them
   into the list at the appropriate positions within the residues.
 
   Usage:
   ======
   hadd(fp,pdb)
   Input:         FILE   *fp        File containing proton generation
                                    parameters.
   Input/Output:  PDB    *pdb       Linked list of protein structure.
   Returns:       int               Number of hydrogens added.
   Externs:       int    info_level Message level.
 
   Revision History:
   =================
   V2.0 16.05.90
   AddH is changed to insert each set of atoms for each PGP, on the
   fly, rather than building a complete list of hydrogens and then
   merging the two lists. This allows us to get round the problem
   of missing atoms, since there will be no merging error.
 
   V2.1 24.05.90
   Returns the number fo hydrogens added. Also fixes bug relating to number
   of type 2 and type 3 H's added.
   Doesn't work under UNIX!
 
   V2.2 15.07.91
   A few bits of tidying up:
   >  Now uses macros.h rather than defining macros itself.
   >  Arrays now changed so should fix alignment problems under UNIX.
   >  Now correctly checks return from forscanf() and no longer
      reads characters to check EOF itself.
   >  Improves treatment of NTER residues where it now generates
      H's on the first true residue. Residues labelled NTER will
      be ignored. Calling FixNterH() will move the H coords into
      the NTER residue if required.
   >  Fixes reading of PGP files with blank lines
   >  Bug fix for type 3's
   Currently untested under UNIX.
 
   V2.3 07.08.91
   Now sets coords of H to 9999.0 if any of the antecedant atoms
   have these coords.

  OMUPD JAW 11/08/92 Changed prototype headers for CR318 
 
***************************************************************************/
 
 
 
/* Global Arrays, etc. */
 
char   gres[MAXTYPE][8], gatom[MAXTYPE][8][8];
char   grname[8], ghname[8][8], gnat[MAXATINRES][8];
int    igtype[MAXTYPE], npgp, no, kmax, it1, it2, it3, it4, it5, ih;
float  gr[MAXTYPE], alpha[MAXTYPE], beta[MAXTYPE];
float  gx[MAXATINRES], gy[MAXATINRES], gz[MAXATINRES], hfac, fac;
 
extern int info_level;
 
/* This is the main entry point */
 
int hadd(
FILE *fp,
PDB  *pdb
)
{
   PDB *p;
   int atomcount=1;
   int nhydrogens;
   unsigned short err_flag=0;
 
   /* Read the parameter file */
   ReadPGP(fp);
 
   /* Generate the hydrogens */
   if((nhydrogens=GenH(pdb,&err_flag,grname))==0)
   {
      prdie("GenH returned an error.\n");
   }
   for(p=pdb;p;NEXT(p)) p->atnum=atomcount++;
   return(nhydrogens);
}
 

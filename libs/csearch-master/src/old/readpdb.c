#include "ProtoTypes.h"
#include "CongenProto.h"
 
/***************************************************************************
 
   Program:
   File:       ReadPDB.c
 
   Version:    V1.2.3
   Date:       04.11.88
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
 
   Written:    While at Laboratory of Molecular Biophysics,
                        University of Oxford,
                        The Rex Richards Building,
                        South Parks Road,
                        Oxford,
                        OX1 3QU.
 
****************************************************************************
 
   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission from
   the author, although it may be given away free with commercial products,
   providing it is made clear that this program is free and that the source
   code is provided with the program.
 
****************************************************************************
 
   Description:
   ============
   ReadPDB(fp,pdb,natom) - This subroutine will read a .PDB file
   of any size and form a linked list of the protein structure.
   This list is contained in a linked set of structures of type
   pdb_entry. The strucure is set up by including the file
   "pdb.h". For details of the structure, see this file.
 
   To free the space created by this routine, call unlink_pdb(pdb).
 
   The parameters passed to the subroutine are:
   fp    - A pointer to type FILE in which the .PDB file is stored.
   pdb   - A pointer to type PDB.
   natom - A pointer to type integer in which the number of atoms
           found is stored.
 
   To define a structure list in which to store the protein, the
   user need only include the file "pdb.h", declare a pointer to a
   structure of type PDB using the statement:
       PDB *mypdb;
   and allocate space for the structure using the macro:
       init_pdb(mypdb);
 
NOTE:  Although some of the fields are represented by a single character,
       they are still stored in character arrays.
 
BUGS:  The subroutine cannot read files with Fortran carriage control!
       It just sits there and page faults like crazy.
 
****************************************************************************
 
   Usage:
   ======
   ReadPDB(fp,pdb,natom) - This subroutine will read a .PDB file
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
            PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list
   Output:  int      *natom   Number of atoms read.
 
 
****************************************************************************
 
   Revision History:
   =================
      V1.1     07.02.89
      Now ignores any records from the .PDB file which don't start
      with ATOM or HETATM.
 
      V1.2     28.03.90
      Some fields altered to match the exact specifications of the PDB.
      The only differences from the standard are:
      1. The residue name is 4 characters rather than 3 (allowing LYSH,
      HISA, etc.).
      2. The atom name starts one column later than the standard and is
      four columns wide encompasing the standard's `alternate' field.
         These two differences from the standard reflect the common
      usage.
 
      V1.2.1   28.06.90
      Buffer size increased to 85 chars.
 
      V1.2.2   15.02.91
      Simply changed comment header to match new standard.
 
      V1.2.3   24.07.91
      Some characters were being set to NULL, now set to '\0'
***************************************************************************/
 
 
#define CR 13
#define LF 10
 
void ReadPDB(
FILE *fp,
PDB *pdb,
int *natom
)
{
   int c,i;
   int SIZE=sizeof(PDB);
   char *buffer,*start;
   char temp[15];
   PDB *p,*last;
 
   *natom = 0;
 
   start = (char *)malloc(85);
 
   p = pdb;
   c = getc(fp);
   while(c != EOF)
   {
      ungetc(c,fp);
 
      /* Set the value for buffer */
      buffer = start;
 
      /* Read a line from the file into buffer[] */
      i=0;
      while((c != LF) && (c != CR))
      {
         c = getc(fp);
         buffer[i] = (char) c;
         i++;
      }
 
      if((!strncmp(buffer,"ATOM  ",6)) || (!strncmp(buffer,"HETATM",6)))
      {
 
         /* Increment the number of atoms */
         (*natom)++;
 
         /* Copy the first 6 charcters into JUNK */
         strncpy(p->junk,buffer,6);
         p->junk[6] = '\0';
 
         buffer += 6;
 
         /* Copy the next 5 characters into temp & convert to integer */
         strncpy(temp,buffer,5);
         temp[5] = '\0';
         p->atnum = atoi(temp);
 
         buffer += 7; /* 2 spaces */
 
         /* Copy the next 4 characters into ATNAM */
         strncpy(p->atnam,buffer,4);
         p->atnam[4] = '\0';
 
         buffer += 4;
 
         /* Copy the next 4 characters into RESNAM */
         strncpy(p->resnam,buffer,4);
         p->resnam[4] = '\0';
 
         buffer += 4;
 
         /* Copy the next 1 character into CHAIN */
         strncpy(p->chain,buffer,1);
         p->chain[1] = '\0';
 
         buffer += 1;
 
         /* Copy the next 4 characters into temp and convert to integer RESNUM */
         strncpy(temp,buffer,4);
         temp[4]='\0';
         p->resnum = atoi(temp);
 
         buffer += 4;
 
 
         /* Copy the next character into INSERT */
         strncpy(p->insert,buffer,1);
         p->insert[1] = '\0';
 
         buffer += 4;
 
         /* Copy the next 8 characters into temp and convert to real X */
         strncpy(temp,buffer,8);
         temp[8]='\0';
         p->x = atof(temp);
 
         buffer += 8;
 
         /* Copy the next 8 characters into temp and convert to real Y */
         strncpy(temp,buffer,8);
         temp[8]='\0';
         p->y = atof(temp);
 
         buffer += 8;
 
         /* Copy the next 8 characters into temp and convert to real Z */
         strncpy(temp,buffer,8);
         temp[8]='\0';
         p->z = atof(temp);
 
         buffer += 8;
 
         /* Check that OCC and BVAL have been specified */
         if(buffer-start < i)
         {
 
            /* Copy the next 6 characters into temp and convert to real OCC */
            strncpy(temp,buffer,6);
            temp[6]='\0';
            p->occ = atof(temp);
 
            buffer += 6;
 
            if(buffer-start < i)
            {
               /* Copy the next 6 characters into temp and convert to real BVAL */
               strncpy(temp,buffer,6);
               temp[6]='\0';
               p->bval = atof(temp);
 
            }
            else
            {
               p->bval=0.0;
            }
         }
         else
         {
            p->occ=0.0;
            p->bval=0.0;
         }
 
         /* Allocate space for the next record */
         p->next = (PDB *)malloc(SIZE);
 
 
         /* and move the pointer on to the next record */
         last = p;
         p = p->next;
 
 
         /* Check we have managed to free sufficient space */
         if (p == NULL)
         {
            printf("Insufficient free memory\n");
            exit(1);
         }
      }
 
      c = getc(fp);
   }
 
   /* Free the extra space made for the non-existant extra record */
   free((char *)p);
   /* and set the ->next pointer from the last record to null */
   last->next = NULL;
}
 

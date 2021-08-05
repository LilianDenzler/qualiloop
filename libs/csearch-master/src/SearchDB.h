#ifndef __SEARCH_DB_H__
#define __SEARCH_DB_H__
/*****************************************************************************
 *      Name:                                                                *
 *  Function:                                                                *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author:                                                                *
 *      Date:                                                                *
 *---------------------------------------------------------------------------*
 *    Inputs:                                                                *
 *   Outputs:                                                                *
 *   Returns:                                                                *
 * Externals:                                                                *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/

#define MAXCONS      20  /* The maximum number of loop constraints */
#define MAXLOOPLEN   20  /* The maximum number of loop residues */
#define MAX_LINE_LEN 80  /* The maximum length of an input line */
#define MAX_FILE_LEN 80  /* The maximum length of a file name */
#define MAX_DIR_LEN  80  /* The maximum length of a directory name */
#define MAXAMINO     20  /* The maximum number of different amino acids */
#define MAXRESAT     17  /* The maximum number of atoms per residue */
#define MAXLOOPAT    340 /* The maximum number of atoms per loop */
#define MAXID        10  /* The maximum loop name character string length */
#define MAXSEQ       25  /* The maximum loop sequence string length */

/* NOTE: The above value of MAXLOOPAT should be MAXRESAT*MAXLOOPLEN */

#define RAD_TO_DEG (float)57.29577951 /* Radian to degrees conversion factor */
#define ZERO       (float)0.0
#define ONE        (float)1.0
#define TWO        (float)2.0
#define TEN        (float)10.0

#define POSITIVE "DP"    /* Identifier for a positive (forward) constraint */
#define NEGATIVE "DM"    /* Identifier for a negative (backward) constraint */
#define COMMENT  "!"     /* A file comment line character */

typedef struct CadEntry CadEntry; /* Structure for C-alpha database entries */
struct CadEntry
{
   CadEntry *next;                    /* Pointer to the next entry */
   CadEntry *last;                    /* Pointer to the previous entry */
   long      Nentry;                  /* The C-alpha database entry number */
   char      PDBcode[5];              /* The entry PDB code */
   char      Chain[4];                /* The entry chain identifier */
   long      ResNum;                  /* The entry residue number */
   char      Insert[4];               /* The entry insert code */
   char      AAcode[2];               /* The entry amino acid code */
   float     Phi;                     /* The amino acid phi torsion angle */
   float     Psi;                     /* The amino acid psi torsion angle */
   float     Omega;                   /* The amino acid omega torsion angle */
   float     CDisPrev[MAXLOOPLEN];    /* The distances to previous C-alphas */
   float     CDisNext[MAXLOOPLEN];    /* The distances to the next C-alphas */
};

typedef struct CadHit CadHit;  /* Structure containing database hits */
struct CadHit
{
   CadHit *next;
   CadHit *last;
   char    PDBcode[5];                 /* The Hit PDB file code */
   char    Chain[4];                   /* The chain ID */
   long    ResNum;                     /* The hit first residue number */
   char    Insert[4];                  /* The hit insert code */
   char    Seq[MAXSEQ];                /* The hit sequence */
   float   Phis[MAXLOOPLEN-1];         /* The phi torsion angles */
   float   Psis[MAXLOOPLEN-1];         /* The psi torsion angles */
   float   Omegas[MAXLOOPLEN-1];       /* The omega torsion angles */
   float   CAdist[MAXLOOPLEN-1];       /* C-alpha distances from 1st C-alpha */
   short   Qclust;                     /* A cluster flag */
};

typedef struct RESTYPE RESTYPE;    /* Structure containing amino acid info */
struct RESTYPE
{
   char  ResNam[MAX_LINE_LEN+1];           /* The amino acid 3 letter code */
   char  AAcode[MAX_LINE_LEN+1];           /* The amino acid 1 letter code */
   int   Natoms;                           /* The number of atoms */
   float X[MAXRESAT];                      /* The atomic X coordinates */
   float Y[MAXRESAT];                      /* The atomic Y coordinates */
   float Z[MAXRESAT];                      /* The atomic Z coordinates */
   char  AtmNam[MAXRESAT][MAX_LINE_LEN+1]; /* The atom names */
};

typedef struct LOOP LOOP;         /* The generic loop coordinate info */
struct LOOP
{
   int   Natoms;                     /* The number of atoms in the loop */
   int   Nres;                       /* The number of residues in the loop */
   int   Natom[MAXAMINO];            /* The loop 'N' atoms */
   int   Hatom[MAXAMINO];            /* The loop 'H' atoms */
   int   Calpha[MAXAMINO];           /* The loop Calpha atoms */
   int   Catom[MAXAMINO];            /* The loop 'C' atoms */
   int   Oatom[MAXAMINO];            /* The loop 'O' atoms */
   float X[MAXLOOPAT];               /* The loop X atom coordinates */
   float Y[MAXLOOPAT];               /* The loop Y atom coordinates */
   float Z[MAXLOOPAT];               /* The loop Z atom coordinates */
   char  AtmNam[MAXLOOPAT][5];       /* The atom names */
   int   ResNum[MAXLOOPAT];          /* The residue number of each atom */
   int   Qinclude[MAXLOOPAT];        /* Is this atom included in the CG file */
};
 
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#endif

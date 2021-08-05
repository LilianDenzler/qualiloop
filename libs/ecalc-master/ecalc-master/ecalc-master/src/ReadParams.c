/*************************************************************************

   Program:    ECalc
   File:       ReadParams.c
   
   Version:    V1.5.2
   Date:       05.03.21
   Function:   Read a parameters file into global arrays
   
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
   V0.1   25.08.94 Original
   V0.2   22.09.94 Added code to get file from environment variable
   V1.0   30.09.94 First release version
   V1.1   11.10.94 Added residue pseudo-energy parameter code
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Fixed for new version of GetWord()
   V1.5.1 07.01.21 Removed unused variables
   V1.5.2 05.03.21 Skipped

*************************************************************************/
/* Includes
*/
#define READPARAMS_MAIN

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/parse.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

#include "ecalc.h"

/************************************************************************/
/* Defines and macros
*/
#define KEY_BOND     0
#define KEY_ANGLE    1
#define KEY_TORSION  2
#define KEY_IMPROPER 3
#define KEY_NONBOND  4
#define KEY_HBOND    5
#define KEY_RESIDUE  6
#define NKEYS        7
#define MAXSTRPARAM  10
#define MAXREALPARAM 1
#define MAXSTRLEN    16

#define CONST        ((REAL)362.3461)

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static BOOL StoreBondParam(char *atom1, char *atom2, REAL p1, REAL p2);
static BOOL StoreAngleParam(char *atom1, char *atom2, char *atom3,
                       REAL p1, REAL p2);
static BOOL StoreTorsionParam(char *atom1, char *atom2, char *atom3, 
                         char *atom4, REAL p1, REAL p2, REAL p3);
static BOOL StoreImproperParam(char *atom1, char *atom2, char *atom3, 
                          char *atom4, REAL p1, REAL p2);
static BOOL StoreNonBondParam(char *atom, int key, int mass, REAL pol, 
                         REAL NEff, REAL vdwr);
static BOOL StoreHBondParam(char *atom1, char *atom2, REAL p1, REAL p2);
static BOOL CalcNBondParams(void);
static BOOL StoreResidueParam(char restype, REAL maxval, REAL scale, 
                              char *distal, int nbox, char *buffer);

/************************************************************************/
/*>BOOL ReadParams(char *filename)
   -------------------------------
   Read a parameter file specified by name

   25.08.94 Original    By: ACRM
   22.09.94 Changed to call OpenFile() rather than fopen()
   11.10.94 Added RESIDUE keyword handling
*/
BOOL ReadParams(char *filename)
{
   KeyWd keywords[NKEYS];
   char  *StrParam[MAXSTRPARAM],
         buffer[MAXBUFF];
   REAL  RealParam[MAXREALPARAM];
   int   i;
   FILE  *fp;
   BOOL  NoVar;
   
   /* Open the file                                                     */
   if((fp=OpenFile(filename,DATADIR,"r",&NoVar)) == NULL)
   {
      sprintf(gError,"Unable to open file `%s'",filename);
      StoreError("ReadParams()",gError);

      if(NoVar)
      {
         sprintf(gError,"Environment variable `%s' not set",DATADIR);
         StoreError("ReadParams()",gError);
      }
      
      return(FALSE);
   }

   /* Allocate memory for return strings                                */
   for(i=0; i<MAXSTRPARAM; i++)
   {
      if((StrParam[i] = (char *)malloc(MAXSTRLEN * sizeof(char))) == NULL)
      {
         StoreError("ReadParams()","No memory for parser");
         return(FALSE);
      }
   }
   
   /* Build keywords                                                    */
   MAKEKEY(keywords[KEY_BOND],     "BOND",     STRING, 4);
   MAKEKEY(keywords[KEY_ANGLE],    "ANGLE",    STRING, 5);
   MAKEKEY(keywords[KEY_TORSION],  "TORSION",  STRING, 7);
   MAKEKEY(keywords[KEY_IMPROPER], "IMPROPER", STRING, 6);
   MAKEKEY(keywords[KEY_NONBOND],  "NONBOND",  STRING, 6);
   MAKEKEY(keywords[KEY_HBOND],    "HBOND",    STRING, 4);
   MAKEKEY(keywords[KEY_RESIDUE],  "RESIDUE",  STRING, 5);

   /* Parse the file                                                    */
   while(fgets(buffer,MAXBUFF,fp))
   {
      switch(parse(buffer,NKEYS,keywords,RealParam,StrParam))
      {
      case PARSE_ERRC:
         sprintf(gError,"Unrecognised keyword `%s'\n",buffer);
         StoreError("ReadParams()",gError);
         break;
      case PARSE_ERRP:
         sprintf(gError,"Error in parameters `%s'\n",buffer);
         StoreError("ReadParams()",gError);
         break;
      case KEY_BOND:
         if(!StoreBondParam(StrParam[0], StrParam[1], 
                            (REAL)atof(StrParam[2]),
                            (REAL)atof(StrParam[3])))
            return(FALSE);
         break;
      case KEY_ANGLE:
         if(!StoreAngleParam(StrParam[0], StrParam[1], StrParam[2],
                             (REAL)atof(StrParam[3]),
                             (REAL)atof(StrParam[4])))
            return(FALSE);
         break;
      case KEY_TORSION:
         if(!StoreTorsionParam(StrParam[0], StrParam[1], StrParam[2], 
                               StrParam[3],
                               (REAL)atof(StrParam[4]),
                               (REAL)atof(StrParam[5]),
                               (REAL)atof(StrParam[6])))
            return(FALSE);
         break;
      case KEY_IMPROPER:
         if(!StoreImproperParam(StrParam[0], StrParam[1], StrParam[2], 
                                StrParam[3],
                                (REAL)atof(StrParam[4]),
                                (REAL)atof(StrParam[5])))
            return(FALSE);
         break;
      case KEY_NONBOND:
         if(!StoreNonBondParam(StrParam[0], atoi(StrParam[1]),
                               (int)atoi(StrParam[2]),
                               (REAL)atof(StrParam[3]),
                               (REAL)atof(StrParam[4]),
                               (REAL)atof(StrParam[5])))
            return(FALSE);
         break;
      case KEY_HBOND:
         if(!StoreHBondParam(StrParam[0], StrParam[1],
                             (REAL)atof(StrParam[2]),
                             (REAL)atof(StrParam[3])))
            return(FALSE);
         break;
      case KEY_RESIDUE:
         if(!StoreResidueParam(StrParam[0][0], 
                               (REAL)atof(StrParam[1]),
                               (REAL)atof(StrParam[2]),
                               StrParam[3],
                               (int)atoi(StrParam[4]),
                               buffer))
            return(FALSE);
         break;
      }
   }
   
   /* Free memory for return strings                                    */
   for(i=0; i<MAXSTRPARAM; i++)
      free(StrParam[i]);

   return(CalcNBondParams());
}
   
/************************************************************************/
/*>static BOOL StoreBondParam(char *atom1, char *atom2, REAL force, 
                              REAL OptLen)
   ----------------------------------------------------------------

   25.08.94 Original    By: ACRM
*/
static BOOL StoreBondParam(char *atom1, char *atom2, REAL force, 
                           REAL OptLen)
{
   if(gNumBondParams >= MAXBONDP)
   {
      StoreError("StoreBondParam()","Too many bond parameters");
      return(FALSE);
   }
   
   strcpy(gBondParams[gNumBondParams].atom1, atom1);
   padterm(gBondParams[gNumBondParams].atom1, 4);
   strcpy(gBondParams[gNumBondParams].atom2, atom2);
   padterm(gBondParams[gNumBondParams].atom2, 4);

   gBondParams[gNumBondParams].force = force;
   gBondParams[gNumBondParams].OptLen = OptLen;

   gNumBondParams++;

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreAngleParam(char *atom1, char *atom2, char *atom3,
                               REAL force, REAL OptLen)
   -------------------------------------------------------------

   25.08.94 Original    By: ACRM
*/
static BOOL StoreAngleParam(char *atom1, char *atom2, char *atom3,
                            REAL force, REAL OptAng)
{
   if(gNumAngleParams >= MAXANGLEP)
   {
      StoreError("StoreAngleParam()","Too many angle parameters");
      return(FALSE);
   }
   
   strcpy(gAngleParams[gNumAngleParams].atom1, atom1);
   padterm(gAngleParams[gNumAngleParams].atom1, 4);
   strcpy(gAngleParams[gNumAngleParams].atom2, atom2);
   padterm(gAngleParams[gNumAngleParams].atom2, 4);
   strcpy(gAngleParams[gNumAngleParams].atom3, atom3);
   padterm(gAngleParams[gNumAngleParams].atom3, 4);

   gAngleParams[gNumAngleParams].force  = force;
   gAngleParams[gNumAngleParams].OptAng = OptAng * PI / (REAL)180.0;

   gNumAngleParams++;

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreTorsionParam(char *atom1, char *atom2, 
                                 char *atom3, char *atom4, REAL force, 
                                 REAL period, REAL OptTor)
   -------------------------------------------------------------------

   25.08.94 Original    By: ACRM
   02.09.94 Added checking of periodicity values and opt tor
*/
static BOOL StoreTorsionParam(char *atom1, char *atom2, char *atom3, 
                              char *atom4, REAL force, REAL period, 
                              REAL OptTor)
{
   int  IPeriod;
   REAL OptTorRad;
   
   
   if(gNumTorsionParams >= MAXTORSIONP)
   {
      StoreError("StoreTorsionParam()","Too many torsion parameters");
      return(FALSE);
   }
   
   IPeriod   = (int)(period+0.0001);
   if(IPeriod==(int)(period-0.0001) || 
      IPeriod>6                     || 
      IPeriod<1                     || 
      IPeriod==5)
   {
      sprintf(gError,"Invalid torsion angle periodicity (%f) in \
parameter file",period);
      
      StoreError("StoreTorsionParam()",gError);
      return(FALSE);
   }

   OptTorRad = OptTor * PI / (REAL)180.0;
   if((ABS(OptTorRad)      > (REAL)0.01) &&
      (ABS(PI - OptTorRad) > (REAL)0.01))
   {
      sprintf(gError,"Invalid optimum torsion angle (%f) in parameter \
file",OptTor);
      
      StoreError("StoreTorsionParam()",gError);
      return(FALSE);
   }

   strcpy(gTorsionParams[gNumTorsionParams].atom1, atom1);
   padterm(gTorsionParams[gNumTorsionParams].atom1, 4);
   strcpy(gTorsionParams[gNumTorsionParams].atom2, atom2);
   padterm(gTorsionParams[gNumTorsionParams].atom2, 4);
   strcpy(gTorsionParams[gNumTorsionParams].atom3, atom3);
   padterm(gTorsionParams[gNumTorsionParams].atom3, 4);
   strcpy(gTorsionParams[gNumTorsionParams].atom4, atom4);
   padterm(gTorsionParams[gNumTorsionParams].atom4, 4);

   gTorsionParams[gNumTorsionParams].force  = force;
   gTorsionParams[gNumTorsionParams].period = IPeriod;
   gTorsionParams[gNumTorsionParams].OptTor = OptTorRad;

   gNumTorsionParams++;

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreImproperParam(char *atom1, char *atom2, char *atom3, 
                                  char *atom4, REAL force, REAL OptTor)
   ---------------------------------------------------------------------

   25.08.94 Original    By: ACRM
*/
static BOOL StoreImproperParam(char *atom1, char *atom2, char *atom3, 
                               char *atom4, REAL force, REAL OptTor)
{
   if(gNumImproperParams >= MAXIMPROPERP)
   {
      StoreError("StoreImproperParam()","Too many Improper parameters");
      return(FALSE);
   }
   
   strcpy(gImproperParams[gNumImproperParams].atom1, atom1);
   padterm(gImproperParams[gNumImproperParams].atom1, 4);
   strcpy(gImproperParams[gNumImproperParams].atom2, atom2);
   padterm(gImproperParams[gNumImproperParams].atom2, 4);
   strcpy(gImproperParams[gNumImproperParams].atom3, atom3);
   padterm(gImproperParams[gNumImproperParams].atom3, 4);
   strcpy(gImproperParams[gNumImproperParams].atom4, atom4);
   padterm(gImproperParams[gNumImproperParams].atom4, 4);

   gImproperParams[gNumImproperParams].force  = force;
   gImproperParams[gNumImproperParams].OptTor = OptTor * PI / (REAL)180.0;

   gNumImproperParams++;

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreNonBondParam(char *atom, REAL pol, REAL NEff, 
                                 REAL vdwr)
   --------------------------------------------------------------

   25.08.94 Original    By: ACRM
*/
static BOOL StoreNonBondParam(char *atom, int key, int mass, REAL pol, 
                              REAL NEff, REAL vdwr)
{
   if(gNumAtomParams >= MAXATOMP)
   {
      StoreError("StoreNonBondParam()","Too many NonBond parameters");
      return(FALSE);
   }
   
   strcpy(gAtomParams[gNumAtomParams].atom, atom);
   padterm(gAtomParams[gNumAtomParams].atom, 4);

   gAtomParams[gNumAtomParams].key  = key;
   gAtomParams[gNumAtomParams].mass = mass;
   gAtomParams[gNumAtomParams].pol  = pol;
   gAtomParams[gNumAtomParams].NEff = NEff;
   gAtomParams[gNumAtomParams].vdwr = vdwr;

   gNumAtomParams++;

   return(TRUE);
}

/************************************************************************/
/*>static BOOL StoreHBondParam(char *atom1, char *atom2, REAL EMin, 
                               REAL RMin)
   ----------------------------------------------------------------

   25.08.94 Original    By: ACRM
   13.09.94 Initialises each item in HBond array
*/
static BOOL StoreHBondParam(char *atom1, char *atom2, REAL EMin, 
                            REAL RMin)
{
   int  i,
        j,
        Key1,
        Key2;
   REAL a,
        b;
   static BOOL FirstCall = TRUE;
   char        Atom1[8],
               Atom2[8];

   strcpy(Atom1,atom1);
   padterm(Atom1,4);
   strcpy(Atom2,atom2);
   padterm(Atom2,4);
   
   if(FirstCall)
   {
      if((gHBondParams = (HBONDPARAM **)
          malloc(gNumAtomParams * sizeof(HBONDPARAM *)))==NULL)
      {
         StoreError("StoreHBondParam()","No memory for HBond parameters");
         return(FALSE);
      }
      for(i=0;i<gNumAtomParams;i++)
      {
         if((gHBondParams[i] = (HBONDPARAM *)
             malloc(gNumAtomParams * sizeof(HBONDPARAM)))==NULL)
         {
            for(j=0; j<i; j++)
               free(gHBondParams[j]);
            StoreError("StoreHBondParam()",
                       "No memory for bond parameters");
            return(FALSE);
         }
         for(j=0; j<gNumAtomParams; j++)
         {
            gHBondParams[i][j].AtomD[0] = '\0';
            gHBondParams[i][j].AtomA[0] = '\0';
            gHBondParams[i][j].r12      = (REAL)0.0;
            gHBondParams[i][j].r10      = (REAL)0.0;
         }
      }
      FirstCall = FALSE;
   }
   
   
   if(gNumHBondParams >= (gNumAtomParams * gNumAtomParams))
   {
      StoreError("StoreHBondParam()","Too many HBond parameters");
      return(FALSE);
   }

   Key1 = (-1);
   for(i = 0; i<gNumAtomParams; i++)
   {
      if(!strncmp(gAtomParams[i].atom,Atom1,4))
      {
         Key1 = gAtomParams[i].key;
         break;
      }
   }

   Key2 = (-1);
   for(i = 0; i<gNumAtomParams; i++)
   {
      if(!strncmp(gAtomParams[i].atom,Atom2,4))
      {
         Key2 = gAtomParams[i].key;
         break;
      }
   }

   if(Key1==(-1) || Key2==(-1))
   {
      StoreError("StoreHBondParam()","HBond atom not found in atom type \
list (NONBONDs must appear before\nHBONDs in Params file)");
      return(FALSE);
   }
   
   strcpy(gHBondParams[Key1][Key2].AtomD, Atom1);
   strcpy(gHBondParams[Key1][Key2].AtomA, Atom2);

   b = (EMin/(-(REAL)pow((double)5.0, (double)5.0) /
               (REAL)pow((double)6.0, (double)6.0))) *
       (REAL)pow((double)(RMin*RMin/(REAL)1.2),(double)5.0);
   a = b * RMin * RMin / (REAL)1.2;

   gHBondParams[Key1][Key2].r12 = a;
   gHBondParams[Key1][Key2].r10 = b;

   gNumHBondParams++;

   return(TRUE);
}

/************************************************************************/
/*>static BOOL CalcNBondParams(void)
   ---------------------------------
   Convert the non-bond atom parameters to r6/r12 constants for each atom
   pair.

   25.08.94 Original    By: ACRM
   05.09.94 Modified such that gNonBondParams is a 2D array
*/
static BOOL CalcNBondParams(void)
{
   int  AtomI,
        AtomJ,
/*
        KeyI,
        KeyJ,
*/
        i;
   REAL Numer, Denom,
        a,     b,
        vdw,   r6,
        emin;

   /* Allocate memory for the gNonBondParams array                      */
   if((gNonBondParams = (NONBONDPARAM **)
       malloc(gNumAtomParams * sizeof(NONBONDPARAM *)))==NULL)
   {
      StoreError("CalcNBondParams()",
                 "No memory for Non-bond parameter array");
      return(FALSE);
   }
   for(i=0; i<gNumAtomParams; i++)
   {
      if((gNonBondParams[i] = (NONBONDPARAM *)
          malloc(gNumAtomParams * sizeof(NONBONDPARAM)))==NULL)
      {
         int j;
         for(j=0; j<i; j++)
            free(gNonBondParams[j]);
            
         StoreError("CalcNBondParams()",
                    "No memory for Non-bond parameter array");
         return(FALSE);
      }
   }
   
   for(AtomI=0; AtomI<gNumAtomParams; AtomI++)
   {
      for(AtomJ=0; AtomJ<gNumAtomParams; AtomJ++)
      {
         if(gNumNonBondParams >= MAXNONBONDP)
         {
            StoreError("CalcNBondParams()",
                       "Non-bond Params size exceeded\n");
            return(FALSE);
         }

         vdw = gAtomParams[AtomI].vdwr + gAtomParams[AtomJ].vdwr;
         r6  = (REAL)pow((double)vdw, (double)6.0);
         if(r6 < (REAL)0.000001) r6 = (REAL)0.000001;
         Numer = CONST * gAtomParams[AtomI].pol
                       * gAtomParams[AtomJ].pol;
         if((gAtomParams[AtomI].NEff >= (REAL)0.5) &&
            (gAtomParams[AtomJ].NEff >= (REAL)0.5) &&
            Numer != (REAL)0.0)
         {
            Denom = (REAL)sqrt((double)(gAtomParams[AtomI].pol /
                                        gAtomParams[AtomI].NEff)) +
                    (REAL)sqrt((double)(gAtomParams[AtomJ].pol /
                                        gAtomParams[AtomJ].NEff));
            b     = Numer/Denom;
            emin  = (REAL)(-0.5)*b/r6;
         }
         else if((gAtomParams[AtomI].NEff <= (REAL)0.0) &&
                 (gAtomParams[AtomJ].NEff <= (REAL)0.0))
         {
            emin = (REAL)(-sqrt((double)(gAtomParams[AtomI].NEff *
                                         gAtomParams[AtomJ].NEff)));
            b    = (REAL)(-2.0)*emin*r6;
         }
         else if((gAtomParams[AtomI].NEff <= (REAL)0.0) &&
                 (gAtomParams[AtomJ].NEff >= (REAL)0.5))
         {
            Numer = (REAL)(0.25 * CONST * gAtomParams[AtomJ].pol *
                           sqrt((double)(gAtomParams[AtomJ].pol *
                                         gAtomParams[AtomJ].NEff)) / r6);
            emin  = (REAL)(-sqrt((double)(-gAtomParams[AtomI].NEff * 
                                          Numer)));
            b     = (REAL)(-2.0)*emin*r6;
         }
         else if((gAtomParams[AtomJ].NEff <= (REAL)0.0) &&
                 (gAtomParams[AtomI].NEff >= (REAL)0.5))
         {
            Numer = (REAL)(0.25 * CONST * gAtomParams[AtomI].pol *
                           sqrt((double)(gAtomParams[AtomI].pol *
                                         gAtomParams[AtomI].NEff)) / r6);
            emin  = (REAL)(-sqrt((double)(-gAtomParams[AtomJ].NEff * 
                                          Numer)));
            b     = (REAL)(-2.0)*emin*r6;
         }
         else
         {
            b     = (REAL)0.0;
            emin  = (REAL)0.0;
         }
         
         a = (REAL)0.5 * b * r6;

/*         
         KeyI = gAtomParams[AtomI].key;
         KeyJ = gAtomParams[AtomJ].key;
*/
         strcpy(gNonBondParams[AtomI][AtomJ].atom1,
                gAtomParams[AtomI].atom);
         strcpy(gNonBondParams[AtomI][AtomJ].atom2,
                gAtomParams[AtomJ].atom);
         gNonBondParams[AtomI][AtomJ].r6  = b;
         gNonBondParams[AtomI][AtomJ].r12 = a;

         gNumNonBondParams++;
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>static BOOL StoreResidueParam(char restype, REAL maxval, REAL scale, 
                                 char *distal, int nbox, char *buffer)
   --------------------------------------------------------------------
   Store a residue parameter specified using the line:
   RESIDUE X maxval scale distal nbox box.....

   11.10.94 Original    By: ACRM
   06.02.03 Fixed for new version of GetWord()
*/
static BOOL StoreResidueParam(char restype, REAL maxval, REAL scale, 
                              char *distal, int nbox, char *buffer)
{
   static RESPARAM *p = NULL;
   char            *buff,
                   word[MAXBUFF];
   int             i;
   
   /* Allocate space in the gResParam linked list                       */
   if(gResParam == NULL)
   {
      INIT(gResParam, RESPARAM);
      p = gResParam;
   }
   else
   {
      ALLOCNEXT(p, RESPARAM);
   }
   if(p==NULL)
   {
      StoreError("StoreResidueParam()","No memory for linked list");
      return(FALSE);
   }
   
   /* Copy input parameters into structure                              */
   p->restype = restype;
   p->maxval  = maxval;
   p->scale   = scale;
   p->nbox    = nbox;
   strcpy(p->distal, distal);
   padterm(p->distal, 4);

   /* Allocate space in structure for the values                        */
   if((p->values = (REAL *)malloc(nbox * sizeof(REAL))) == NULL)
   {
      StoreError("StoreResidueParam()","No memory for values");
      return(FALSE);
   }

   /* Eat the first 6 words out of the buffer                           */
   for(i=0, buff=buffer; i<6; i++)
      buff = GetWord(buff, word, MAXBUFF);

   for(i=0; i<p->nbox; i++)
      p->values[i] = (REAL)0.0;

   /* Put the values into the structure                                 */
   for(i=0; buff!=NULL && i<p->nbox; i++)
   {
      buff = GetWord(buff, word, MAXBUFF);
      sscanf(word, "%lf", &(p->values[i]));
   }
      
   return(TRUE);
}

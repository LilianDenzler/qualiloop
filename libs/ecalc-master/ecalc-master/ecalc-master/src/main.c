/*************************************************************************

   Program:    ECalc
   File:       main.c
   
   Version:    V1.5.2
   Date:       05.03.21
   Function:   Main program for energy calculation
   
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
   V0.1   09.09.94 Original
   V1.0   30.09.94 First release version
   V1.1   12.10.94 Added code to handle residue part of potential
                   Fixed bug in low energy caching
   V1.2   10.11.94 Bug fixes in energy.c
   V1.3   28.11.94 Bug fixes in energy.c
   V1.4   18.05.95 Shake support
   V1.5   06.02.03 Bug fixes in ReadStructure.c
   V1.5.1 07.01.21 General tidy up
   V1.5.2 05.03.21 Bug fixes in energy.c, shake.c and StoreZone()

*************************************************************************/
#define ECALC_MAIN

/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "bioplib/parse.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

#include "ecalc.h"

/************************************************************************/
/* Defines and macros
*/
#define KEY_PDBFILE          0
#define KEY_CGFILE           1
#define KEY_OUTPUTFILE       2
#define KEY_POTENTIAL        3
#define KEY_BONDS            4
#define KEY_ANGLES           5
#define KEY_TORSIONS         6
#define KEY_IMPROPERS        7
#define KEY_HBONDS           8
#define KEY_VDWA             9
#define KEY_VDWR            10
#define KEY_ELECTROSTATIC   11
#define KEY_END             12
#define KEY_DISPLAY         13
#define KEY_CACHE           14
#define KEY_GRIDCUT         15
#define KEY_NONBONDCUT      16
#define KEY_ETA             17
#define KEY_CONSTDIELECTRIC 18
#define KEY_DISTDIELECTRIC  19
#define KEY_CUTONHB         20
#define KEY_CUTOFFHB        21
#define KEY_CUTONANG        22
#define KEY_CUTOFFANG       23
#define KEY_SHOWPARAMS      24
#define KEY_SHOWRTOP        25
#define KEY_SHOWTIMINGS     26
#define KEY_PARAMFILE       27
#define KEY_RTOPFILE        28
#define KEY_DEBUG           29
#define KEY_DISULPHIDES     30
#define KEY_RUN             31
#define KEY_REGRID          32
#define KEY_ZONE            33
#define KEY_IGNORE          34
#define KEY_RESIDUE         35
#define KEY_RELAX           36
#define KEY_TOL             37
#define NKEYWORDS           38
#define MAXSTRPARAM          8
#define MAXREALPARAM         4
#define MAXSTRLEN          160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
#include "protos.h"

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for energy calculation.

   09.09.94 Original    By: ACRM
   13.09.94 Added Flags.Debug in call to ShowErrors()
            Output to file now possible
            Fixed code for printing missing atoms
   14.09.94 Removed unused variables
   15.09.94 Prints message if DoEnergyCalculation() fails
   16.09.94 Allocates memory for cache and displays contents
*/
int main(int argc, char **argv)
{
   MOLECULE *mol;
   FILE     *pdbfp     = NULL,
            *controlfp = NULL,
            *conffp    = NULL,
            *out       = stdout;
/*   BOOL     AtomError  = FALSE; */
   EPARAMS  EParams;
   FLAGS    Flags;
   int      NumSS      = 0;
   char     PDBFile[MAXBUFF],
            ControlFile[MAXBUFF],
            ParamFile[MAXBUFF],
            TopFile[MAXBUFF];

   SetDefaults(&EParams, &Flags, ParamFile, TopFile);

   if(ParseCmdLine(argc, argv, PDBFile, ControlFile))
   {
      if(OpenFiles(PDBFile, ControlFile, &pdbfp, &controlfp))
      {
         if(ParseControlFile(controlfp, &EParams, ParamFile, TopFile,
                             &Flags, &conffp, &pdbfp, &out))
         {
            if(ReadParams(ParamFile))
            {
               if(Flags.PrintParams)
                  PrintParams(out);
               
               if(ReadRTop(TopFile))
               {
                  if(Flags.PrintRTop)
                     PrintRTop(out);
                  
                  if((mol = ReadStructure(pdbfp, Flags.AutoDisulphide, 
                                          &gAtomError, &NumSS))!=NULL)
                  {
                     if(Flags.AutoDisulphide)
                        fprintf(out, "%d disulphides were found and \
patched\n",NumSS);
                     else
                        fprintf(out, "No search was made for \
disulphides\n");

                     FreeRTop();
                     
                     if(gAtomError)
                     {
                        int chain, i, resnum;
                        
                        fprintf(out, "There were missing atoms!\n");
                        for(chain=0; chain<mol->NChains; chain++)
                        {
                           for(i=0; i<(mol->topol[chain])->NAtoms; i++)
                           {
                              if((mol->topol[chain])->atoms[i]->x > 
                                  9998.0 &&
                                 (mol->topol[chain])->atoms[i]->y > 
                                  9998.0 &&
                                 (mol->topol[chain])->atoms[i]->z > 
                                  9998.0)
                              {
                                 resnum = 
                                    GetResFromAtomNum(mol->topol[chain], 
                                                      i);
                                 
                                 fprintf(out, "Chain: %c Residue: %c %d \
Atom: %s\n", 
                                    (mol->topol[chain])->ChainName,
                                    (mol->topol[chain])->sequence[resnum],
                                    resnum+1,
                                    (mol->topol[chain])->atoms[i]->atom);
                              }
                           }
                        }
                     }

                     /* Allocate memory for lowest energy cache         */
                     if(!BuildCache(Flags.NCache))
                     {
                        StoreError("main()","No memory for cache");
                     }
                     else
                     {
                        if(DoEnergyCalculations(out, mol, EParams, 
                                                Flags, conffp))
                        {
                           fprintf(out,"\n\nSuccessful completion\n");
                           ShowCache(out, Flags.NCache);
                        }
                        else
                        {
                           fprintf(out,"\n\nErrors in energy \
calculation\n");
                        }
                        
                     }
                  }
                  else
                  {
                     FreeRTop();
                  }
               }
            }
         }
      }         
   }
   else
   {
      UsageExit();
   }
   
   ShowErrors(NULL, Flags.Debug);
   
   return(0);
}

/************************************************************************/
/*>void SetDefaults(EPARAMS *EParams, FLAGS *Flags, char *ParamFile, 
                    char *TopFile)
   -----------------------------------------------------------------
   Set up various default values

   09.09.94 Original    By: ACRM
   13.09.94 Added Flags->AutoDisulphide
   15.09.94 Added Flags->Regrid
   16.09.94 Changed CutOnHB from 3.0 and GridCut from 12.0
   11.10.94 Added residue stuff
   19.05.95 Added Relaxation parameters to EParams
   06.02.03 Changed default flags so that the residue energy isn't
            done by default
*/
void SetDefaults(EPARAMS *EParams, FLAGS *Flags, char *ParamFile, 
                 char *TopFile)
{
   EParams->GridCut         = 30.0;
   EParams->NonBondCut      = 8.0;
   EParams->eta             = 50.0;
   EParams->ConstDielectric = FALSE;
   EParams->CutOnHB         = 4.0;
   EParams->CutOffHB        = 5.0;
   EParams->CutOnHBAng      = 90.0 * PI / 180.0;
   EParams->CutOffHBAng     = 90.0 * PI / 180.0;
   EParams->DoRelax         = FALSE;
   EParams->VDWTol          = (REAL)0.5;
   EParams->ShakeTol        = (REAL)0.001;

   /* Should these parts of the energy be displayed                     */
   Flags->PrintParams       = FALSE;
   Flags->PrintRTop         = FALSE;
   Flags->ShowBonds         = FALSE;
   Flags->ShowAngles        = FALSE;
   Flags->ShowTorsions      = FALSE;
   Flags->ShowImpropers     = FALSE;
   Flags->ShowVDWA          = FALSE;
   Flags->ShowVDWR          = FALSE;
   Flags->ShowElect         = FALSE;
   Flags->ShowHBonds        = FALSE;
   Flags->ShowResidue       = FALSE;

   Flags->Timings           = FALSE;
   Flags->Debug             = FALSE;

   Flags->AutoDisulphide    = TRUE;

   Flags->NCache            = 1;
   Flags->Regrid            = 0;

   /* Should these parts of the energy be calculated                    */
   Flags->bonds             = TRUE;
   Flags->angles            = TRUE;
   Flags->torsions          = TRUE;
   Flags->impropers         = TRUE;
   Flags->vdwa              = TRUE;
   Flags->vdwr              = TRUE;
   Flags->hbonds            = TRUE;
   Flags->elect             = TRUE;
   Flags->residue           = FALSE;

   /* Scale factors for these parts of the energy                       */
   Flags->BondScale         = (REAL)1.0;
   Flags->AngleScale        = (REAL)1.0;
   Flags->TorsionScale      = (REAL)1.0;
   Flags->ImproperScale     = (REAL)1.0;
   Flags->VDWAScale         = (REAL)1.0;
   Flags->VDWRScale         = (REAL)1.0;
   Flags->HBondScale        = (REAL)1.0;
   Flags->ElectScale        = (REAL)1.0;
   Flags->ResidueScale      = (REAL)0.0;

   strcpy(ParamFile, PARAMFILE);
   strcpy(TopFile,   TOPFILE);
}
   
/************************************************************************/
/*>BOOL OpenFiles(char *PDBFile, char *ControlFile, FILE **pdbfp, 
                  FILE **controlfp)
   --------------------------------------------------------------
   Open the PDB and control files if specified

   09.09.94 Original    By: ACRM
   13.09.94 Added option for control on stdin
*/
BOOL OpenFiles(char *PDBFile, char *ControlFile, FILE **pdbfp, 
               FILE **controlfp)
{
   /* Assume no files                                                   */
   *controlfp = NULL;
   *pdbfp     = NULL;

   /* If a PDB file is specified, open it                               */
   if(PDBFile[0])
   {
      if(PDBFile[0] == '-')     /* PDB file will be on stdin            */
      {
         *pdbfp     = stdin;
      }
      else if((*pdbfp = fopen(PDBFile,"r")) == NULL)
      {
         sprintf(gError,"Unable to open pdb file: %s",PDBFile);
         StoreError("OpenFiles()",gError);
         return(FALSE);
      }
   }

   /* Open the control file if specified                                */
   if(ControlFile[0])
   {
      if(ControlFile[0] == '-')   /* COntrol will be on stdin           */
      {
         *controlfp = stdin;
      }
      else if((*controlfp = fopen(ControlFile,"r")) == NULL)
      {
         sprintf(gError,"Unable to open control file: %s",ControlFile);
         StoreError("OpenFiles()",gError);
         
         return(FALSE);
      }
   }

   /* If no control file has been specified, the PDB file must be       */
   if(*controlfp == NULL && *pdbfp == NULL)
   {
      StoreError("OpenFiles()", "You must specify a PDB file if there is \
no control file");
      return(FALSE);
   }

   /* Check that the files aren't both on stdin                         */
   if(*controlfp == stdin && *pdbfp == stdin)
   {
      StoreError("OpenFiles()", "You can't specify both the PDB file and \
the control file as `-'");
      return(FALSE);
   }
   
   return(TRUE);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *PDBFile, 
                     char *ControlFile)
   -------------------------------------------------------
   Parse the command line

   09.09.94 Original    By: ACRM
   13.09.94 Added call to UsageExit()
            Added - flag on its own to indicate control on stdin
   19.02.21 Now returns FALSE  if no arguments instead of calling 
            UsageExit()
*/
BOOL ParseCmdLine(int argc, char **argv, char *PDBFile, char *ControlFile)
{
   argc--;
   argv++;

   PDBFile[0]     = '\0';
   ControlFile[0] = '\0';

   if(!argc)
      return(FALSE);
  
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'p':
            argc--;
            argv++;
            strncpy(PDBFile,argv[0],MAXBUFF);
            PDBFile[MAXBUFF-1] = '\0';
            break;
         case '\0':
            strcpy(ControlFile,"-");
            return(TRUE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there is only 1 argument left                    */
         if(argc > 1)
            return(FALSE);
         
         /* Copy the argument to ControlFile                            */
         strcpy(ControlFile, argv[0]);
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void UsageExit(void)
   --------------------
   Prints a usage and copyright message and exits the program

   13.09.94 Original    By: ACRM
   11.10.94 V1.1
   10.11.94 V1.2
   28.11.94 V1.3
   18.05.95 V1.4
   06.02.03 V1.5
   07.01.21 V1.5.1
   05.03.21 V1.5.2 Updated file preparation message
*/
void UsageExit(void)
{
   fprintf(stderr,"\nECalc V1.5.2 (c) 1994-2021, Dr. Andrew C.R. Martin, \
University College London.\n\n");

   fprintf(stderr,"Calculates the energy of a protein structure or of a \
set of conformations\n");
   fprintf(stderr,"for a region of a structure stored in CSearch \
(CHARMM-free CONGEN) format.\n");
   fprintf(stderr,"The protein structure file must have hydrogens and \
CHARMM-style NTER\n");
   fprintf(stderr,"and CTER residues. This can be achieved using the \
following UNIX command:\n\n");

   fprintf(stderr,"   pdbchain xxxx.pdb | pdbhadd -c | pdbcter -c | \
pdbrenum > xxxx.pdh\n\n");

   fprintf(stderr,"Usage: ecalc [-p (xxxx.pdh | -)] [control.dat | -]\n");
   fprintf(stderr,"              -p specifies the PDB file (with \
hydrogens). This must be\n");
   fprintf(stderr,"                 specified if no control file is \
given.\n");
   fprintf(stderr,"              If a file is specified as - standard \
input will be used.\n\n");

   fprintf(stderr,"See documentation for the format of the control file, \
or use the X-windows\n");
   fprintf(stderr,"interface.\n\n");

   exit(0);
}


/************************************************************************/
/*>BOOL ParseControlFile(FILE *controlfp, EPARAMS *EParams, 
                         char *ParamFile, char *TopFile, FLAGS *Flags, 
                         FILE **conffp, FILE **pdbfp, FILE **outfp)
   -------------------------------------------------------------------
   Parse the control file (if opened). Sets up the parser, handles the
   file and frees resources allocated for the parser.

   12.09.94 Original    By: ACRM
   13.09.94 Added RUN keyword
   15.09.94 Added REGIRD
   16.09.94 Added ZONE
   17.09.94 Added IGNORE. 
            Modified disulphide zone to use StoreZone() & changed calls
            to StoreZone()
   19.05.95 Added RELAX and TOLERENCE
*/
BOOL ParseControlFile(FILE *controlfp, EPARAMS *EParams, char *ParamFile,
                      char *TopFile, FLAGS *Flags, FILE **conffp, 
                      FILE **pdbfp, FILE **outfp)
{
   MKeyWd KeyWords[NKEYWORDS];
   char   *StrParam[MAXSTRPARAM];
   REAL   RealParam[MAXREALPARAM];
   int    i,
          key,
          nparams,
          state = STATE_NONE;   /* Used for POTENTIAL and DISPLAY       */
   char   buffer[MAXBUFF];
   ZONE   *z;

   /* Check there is a control file to parse                            */
   if(controlfp == NULL)
      return(TRUE);

   /* Allocate memory for returned strings array                        */
   for(i=0; i<MAXSTRPARAM; i++)
   {
      if((StrParam[i] = (char *)malloc(MAXSTRLEN * sizeof(char))) == NULL)
      {
         StoreError("ParseControlFile()",
                    "No memory for parser returned string array");
         return(FALSE);
      }
   }
   
   /* Build the keywords                                                */
   MAKEMKEY(KeyWords[KEY_PDBFILE],         "PDBFILE",         STRING,1,1);
   MAKEMKEY(KeyWords[KEY_CGFILE],          "CONFFILE",        STRING,1,1);
   MAKEMKEY(KeyWords[KEY_OUTPUTFILE],      "OUTFILE",         STRING,1,1);
   MAKEMKEY(KeyWords[KEY_POTENTIAL],       "POTENTIAL",       STRING,0,0);
   MAKEMKEY(KeyWords[KEY_BONDS],           "BONDS",           NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_ANGLES],          "ANGLES",          NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_TORSIONS],        "TORSIONS",        NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_IMPROPERS],       "IMPROPERS",       NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_HBONDS],          "HBONDS",          NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_VDWA],            "VDWA",            NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_VDWR],            "VDWR",            NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_ELECTROSTATIC],   "ELECTROSTATICS",  NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_END],             "END",             STRING,0,0);
   MAKEMKEY(KeyWords[KEY_DISPLAY],         "DISPLAY",         STRING,0,0);
   MAKEMKEY(KeyWords[KEY_CACHE],           "CACHE",           NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_GRIDCUT],         "GRIDCUT",         NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_NONBONDCUT],      "NONBONDCUT",      NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_ETA],             "ETA",             NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_CONSTDIELECTRIC], "CONSTDIELECTRIC", STRING,0,0);
   MAKEMKEY(KeyWords[KEY_DISTDIELECTRIC],  "DISTDIELECTRIC",  STRING,0,0);
   MAKEMKEY(KeyWords[KEY_CUTONHB],         "CUTONHB",         NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_CUTOFFHB],        "CUTOFFHB",        NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_CUTONANG],        "CUTONANG",        NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_CUTOFFANG],       "CUTOFFANG",       NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_SHOWPARAMS],      "SHOWPARAMS",      STRING,0,0);
   MAKEMKEY(KeyWords[KEY_SHOWRTOP],        "SHOWRTOP",        STRING,0,0);
   MAKEMKEY(KeyWords[KEY_SHOWTIMINGS],     "SHOWTIMINGS",     STRING,0,0);
   MAKEMKEY(KeyWords[KEY_PARAMFILE],       "PARAMFILE",       STRING,1,1);
   MAKEMKEY(KeyWords[KEY_RTOPFILE],        "RTOPFILE",        STRING,1,1);
   MAKEMKEY(KeyWords[KEY_DEBUG],           "DEBUG",           STRING,0,0);
   MAKEMKEY(KeyWords[KEY_DISULPHIDES],     "DISULPHIDES",     STRING,1,2);
   MAKEMKEY(KeyWords[KEY_RUN],             "RUN",             STRING,0,0);
   MAKEMKEY(KeyWords[KEY_REGRID],          "REGRID",          NUMBER,1,1);
   MAKEMKEY(KeyWords[KEY_ZONE],            "ZONE",            STRING,2,2);
   MAKEMKEY(KeyWords[KEY_IGNORE],          "IGNORE",          STRING,2,3);
   MAKEMKEY(KeyWords[KEY_RESIDUE],         "RESIDUE",         NUMBER,0,1);
   MAKEMKEY(KeyWords[KEY_RELAX],           "RELAX",           NUMBER,0,0);
   MAKEMKEY(KeyWords[KEY_TOL],             "TOLERENCE",       STRING,2,2);

   /* Check that the last memory allocation succeeded                   */
   if(KeyWords[KEY_DISULPHIDES].name == NULL)
   {
      StoreError("ParseControlFile()",
                 "No memory for parser KeyWords array");
      return(FALSE);
   }

   /* Sit in main parser loop                                           */
   while(fgets(buffer, MAXBUFF, controlfp))
   {
      TERMINATE(buffer);
      
      key = mparse(buffer,NKEYWORDS,KeyWords,RealParam,StrParam,
                   &nparams);
      if(key == KEY_RUN)
         break;
      
      switch(key)
      {
      case PARSE_ERRC:
         fprintf(stderr, "Unknown command: %s (ignored)\n", buffer);
         break;
      case PARSE_ERRP:
         fprintf(stderr, "Error in parameters: %s (ignored)\n", buffer);
         break;
      case KEY_PDBFILE:
         if((*pdbfp = fopen(StrParam[0],"r")) == NULL)
         {
            sprintf(gError,"Unable to open PDB file %s specified in \
control file",StrParam[0]);
            StoreError("ParseControlFile()", gError);
            return(FALSE);
         }
         break;
      case KEY_CGFILE:
         if((*conffp = fopen(StrParam[0],"r")) == NULL)
         {
            sprintf(gError,"Unable to open conformation file %s \
specified in control file",StrParam[0]);
            StoreError("ParseControlFile()", gError);
            return(FALSE);
         }
         break;
      case KEY_OUTPUTFILE:
         if((*outfp = fopen(StrParam[0],"w")) == NULL)
         {
            sprintf(gError,"Unable to open output file %s specified in \
control file",StrParam[0]);
            StoreError("ParseControlFile()", gError);
            return(FALSE);
         }
         break;
      case KEY_POTENTIAL:
         state = STATE_POTENTIAL;
         Flags->bonds     = FALSE;
         Flags->angles    = FALSE;
         Flags->torsions  = FALSE;
         Flags->impropers = FALSE;
         Flags->hbonds    = FALSE;
         Flags->elect     = FALSE;
         Flags->vdwa      = FALSE;
         Flags->vdwr      = FALSE;
         Flags->residue   = FALSE;
         break;
      case KEY_BONDS:
         CHECK_STATE(KeyWords[KEY_BONDS].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowBonds = TRUE;
         }
         else
         {
            Flags->bonds = TRUE;
            if(nparams)
               Flags->BondScale = RealParam[0];
         }
         break;
      case KEY_ANGLES:
         CHECK_STATE(KeyWords[KEY_ANGLES].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowAngles = TRUE;
         }
         else
         {
            Flags->angles = TRUE;
            if(nparams)
               Flags->AngleScale = RealParam[0];
         }
         break;
      case KEY_TORSIONS:
         CHECK_STATE(KeyWords[KEY_TORSIONS].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowTorsions = TRUE;
         }
         else
         {
            Flags->torsions = TRUE;
            if(nparams)
               Flags->TorsionScale = RealParam[0];
         }
         break;
      case KEY_IMPROPERS:
         CHECK_STATE(KeyWords[KEY_IMPROPERS].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowImpropers = TRUE;
         }
         else
         {
            Flags->impropers = TRUE;
            if(nparams)
               Flags->ImproperScale = RealParam[0];
         }
         break;
      case KEY_HBONDS:
         CHECK_STATE(KeyWords[KEY_HBONDS].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowHBonds = TRUE;
         }
         else
         {
            Flags->hbonds = TRUE;
            if(nparams)
               Flags->HBondScale = RealParam[0];
         }
         break;
      case KEY_RESIDUE:
         CHECK_STATE(KeyWords[KEY_RESIDUE].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowResidue = TRUE;
         }
         else
         {
            Flags->residue = TRUE;
            if(nparams)
               Flags->ResidueScale = RealParam[0];
         }
         break;
      case KEY_VDWA:
         CHECK_STATE(KeyWords[KEY_VDWA].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowVDWA = TRUE;
         }
         else
         {
            Flags->vdwa = TRUE;
            if(nparams)
               Flags->VDWAScale = RealParam[0];
         }
         break;
      case KEY_VDWR:
         CHECK_STATE(KeyWords[KEY_VDWR].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowVDWR = TRUE;
         }
         else
         {
            Flags->vdwr = TRUE;
            if(nparams)
               Flags->VDWRScale = RealParam[0];
         }
         break;
      case KEY_ELECTROSTATIC:
         CHECK_STATE(KeyWords[KEY_ELECTROSTATIC].name);
         if(state == STATE_SHOW)
         {
            Flags->ShowElect = TRUE;
         }
         else
         {
            Flags->elect = TRUE;
            if(nparams)
               Flags->ElectScale = RealParam[0];
         }
         break;
      case KEY_END:
         CHECK_STATE(KeyWords[KEY_END].name);
         state = STATE_NONE;
         break;
      case KEY_DISPLAY:
         state = STATE_SHOW;
         break;
      case KEY_CACHE:
         Flags->NCache = (int)RealParam[0];
         break;
      case KEY_GRIDCUT:
         EParams->GridCut = RealParam[0];
         break;
      case KEY_NONBONDCUT:
         EParams->NonBondCut = RealParam[0];
         break;
      case KEY_ETA:
         EParams->eta = RealParam[0];
         break;
      case KEY_CONSTDIELECTRIC:
         EParams->ConstDielectric = TRUE;
         break;
      case KEY_DISTDIELECTRIC:
         EParams->ConstDielectric = FALSE;
         break;
      case KEY_CUTONHB:
         EParams->CutOnHB = RealParam[0];
         break;
      case KEY_CUTOFFHB:
         EParams->CutOffHB = RealParam[0];
         break;
      case KEY_CUTONANG:
         EParams->CutOnHBAng = (PI * RealParam[0]) / (REAL)180.0;
         break;
      case KEY_CUTOFFANG:
         EParams->CutOffHBAng = (PI * RealParam[0]) / (REAL)180.0;
         break;
      case KEY_SHOWPARAMS:
         Flags->PrintParams = TRUE;
         break;
      case KEY_SHOWRTOP:
         Flags->PrintRTop = TRUE;
         break;
      case KEY_SHOWTIMINGS:
         Flags->Timings = TRUE;
         break;
      case KEY_PARAMFILE:
         strcpy(ParamFile, StrParam[0]);
         break;
      case KEY_RTOPFILE:
         strcpy(TopFile, StrParam[0]);
         break;
      case KEY_DEBUG:
         Flags->Debug = TRUE;
         break;
      case KEY_DISULPHIDES:
         Flags->AutoDisulphide = FALSE;
         if(!upstrncmp(StrParam[0], "OFF", 3))
            break;
         if(nparams == 2)
         {
            if(!StoreZone(&gUserSS, StrParam[0], StrParam[1]))
            {
               sprintf(gError,"DISULPHIDE command (%s) in control file \
failed", buffer);
               StoreError("ParseControlFile()", gError);
               return(FALSE);
            }
         }
         else
         {
            sprintf(gError,"Invalid DISULPHIDE command (%s) in command \
file", buffer);
            StoreError("ParseControlFile()", gError);
            return(FALSE);
         }
         break;
      case KEY_REGRID:
         Flags->Regrid = (int)RealParam[0];
         break;
      case KEY_ZONE:
         if(!StoreZone(&gZone, StrParam[0], StrParam[1]))
         {
            StoreError("ParseControlFile()", "Unable to store zone");
            return(FALSE);
         }
         break;
      case KEY_IGNORE:
         if((z = StoreZone(&gIgnore, StrParam[0], StrParam[1]))==NULL)
         {
            StoreError("ParseControlFile()", 
                       "Unable to store ignore zone");
            return(FALSE);
         }
         
         if(nparams == 3)
         {
            if(toupper(StrParam[2][0]) == 'S')
            {
               z->sc = TRUE;
            }
            else
            {
               sprintf(gError,"Unknown IGNORE option (%s) in control \
file", StrParam[2]);
               StoreError("ParseControlFile()", gError);
               return(FALSE);
            }
         }
         break;
      case KEY_RELAX:
         EParams->DoRelax = TRUE;
         break;
      case KEY_TOL:
         switch(toupper(StrParam[0][0]))
         {
         case 'S':
            EParams->ShakeTol = atof(StrParam[1]);
            break;
         case 'V':
            EParams->VDWTol   = atof(StrParam[1]);
            break;
         default:
            sprintf(gError,"Unknown TOLERENCE type (%s) in control \
file", StrParam[0]);
            StoreError("ParseControlFile()", gError);
            return(FALSE);
         }
         break;
      default:
         break;
      }
   }  /* End of main parser loop                                        */
   
   /* Free memory from the KeyWords                                     */
   for(i=0; i<NKEYWORDS; i++)
      if(KeyWords[i].name != NULL) free(KeyWords[i].name);
      
   /* Free the returned string array                                    */
   for(i=0; i<MAXSTRPARAM; i++)
      free(StrParam[i]);
   
   /* Check that a PDB file was specified somewhere (either open before
      calling this routine, or opened by it)
   */
   if(*pdbfp == NULL)
   {
      StoreError("ParseControlFile()","No PDB file has been specified");
      return(FALSE);
   }

   return(TRUE);
}

/************************************************************************/
/*>ZONE *StoreZone(ZONE **zone, char *resspec1, char *resspec2)
   -----------------------------------------------------------
   Stores a zone over which the energy will be calculated. Specifications
   are added to a global linked list

   16.09.94 Original    By: ACRM
   21.09.94 Added zone as parameter rather than global. 
            Corrected error messages
            Returns a pointer to the latest zone rather than BOOL
   05.03.21 Fixed strncpy length
*/
ZONE *StoreZone(ZONE **zone, char *resspec1, char *resspec2)
{
   ZONE *p;
   char chain1,
        chain2,
        insert1,
        insert2;
   int  resnum1,
        resnum2;
   
   /* Split up the residue specifications                               */
   ParseResSpec(resspec1, &chain1, &resnum1, &insert1);
   ParseResSpec(resspec2, &chain2, &resnum2, &insert2);

   /* Check that consecutive numbering has been used                    */
   if(insert1 != ' ' || insert2 != ' ')
   {
      StoreError("StoreZone()",
                 "Consecutive numbering (no insertions) must be used");
      return(NULL);
   }

   /* Check that both res specs are in the same chain                   */
   if(chain1 != chain2)
   {
      StoreError("StoreZone()",
                 "Both residues must be in the same chain");
      return(NULL);
   }

   if(*zone == NULL)
   {
      INIT((*zone), ZONE);
      p = (*zone);
   }
   else
   {
      p = (*zone);
      LAST(p);
      ALLOCNEXT(p, ZONE);
   }
   
   if(p==NULL)
   {
      StoreError("StoreZone()",
                 "No memory to store zone data");
      return(NULL);
   }

   /* 05.03.21 Fixed from 16                                            */
   strncpy(p->res1, resspec1, 15);
   strncpy(p->res2, resspec2, 15);

   p->res1[15] = '\0';
   p->res2[15] = '\0';
   p->sc       = FALSE;

   return(p);
}


/************************************************************************/
/*>BOOL BuildCache(int size)
   -------------------------
   Allocates memory and zeros the lowest energy conformation cache

   16.09.94 Original    By: ACRM
*/
BOOL BuildCache(int size)
{
   int i;
  
   /* If size is 0 (or less) just return                                */
   if(size < 1)
      return(TRUE);
 
   if((gCache = (CACHE *)malloc(size * sizeof(CACHE)))==NULL)
      return(FALSE);

   for(i=0; i<size; i++)
      gCache[i].confnum = 0;

   return(TRUE);
}

      
/************************************************************************/
/*>void CacheConf(int ConfNum, REAL energy, int CacheSize)
   -------------------------------------------------------
   See if a conformation neds stashing in the low energy cache (gCache)
   and do so if necessary.

   16.09.94 Original    By: ACRM
   19.09.94 Added check that there is a cache (!)
   12.10.94 Sets CacheHigh to zero before re-scanning for highest
*/
void CacheConf(int ConfNum, REAL energy, int CacheSize)
{
   static int  CacheUsed = 0,
               CacheHigh;
   static REAL HighestE;
   int         i;

   /* Just return if there is no cache                                  */
   if(CacheSize < 1)
      return;
   
   /* Just add to the cache until it is full                            */
   if(CacheUsed < CacheSize)
   {
      gCache[CacheUsed].energy  = energy;
      gCache[CacheUsed].confnum = ConfNum;

      if(CacheUsed == 0)
      {
         CacheHigh = 0;
         HighestE  = energy;
      }
      else if(energy > HighestE)
      {
         CacheHigh = CacheUsed;
         HighestE  = energy;
      }

      CacheUsed++;
   }
   else if(energy < HighestE)
   {
      /* Store this conformation in the cache                           */
      gCache[CacheHigh].energy  = energy;
      gCache[CacheHigh].confnum = ConfNum;
      
      /* Refind the highest energy in the cache                         */
      HighestE  = gCache[0].energy;
      CacheHigh = 0;
      for(i=1; i<CacheSize; i++)
      {
         if(gCache[i].energy > HighestE)
         {
            CacheHigh = i;
            HighestE  = gCache[i].energy;
         }
      }
   }
   
   return;
}

/************************************************************************/
/*>void ShowCache(FILE *out, int NCache)
   -------------------------------------
   Displays the contents of the low energy conformation cache

   16.09.94 Original    By: ACRM
*/
void ShowCache(FILE *out, int NCache)
{
   int  i,
        NConfs;
   
   if(NCache)
   {
      for(i=0, NConfs=0; i<NCache; i++)
      {
         if(gCache[i].confnum)
            NConfs++;
      }

      if(NConfs)
      {
         if(NConfs == 1)
            fprintf(out, "\nThe lowest energy conformation was:\n\n");
         else
            fprintf(out, "\nThe lowest energy %d conformations were:\n\n",
                    NConfs);
      }
      
      for(i=0; i<NCache; i++)
      {
         if(gCache[i].confnum)
         {
            fprintf(out, "Conformation %5d : Energy %f\n",
                    gCache[i].confnum, gCache[i].energy);
         }
      }
   }
}
   
      


#ifdef JUNK

/************************************************************************/
/*>BOOL StoreUserDisulphide(char *resspec1, char *resspec2)
   --------------------------------------------------------
   Adds the pair of residue specifications to the global linked list
   of user disulphide specifications

   13.09.94 Original    By: ACRM
   16.09.94 Changed name of USERSS to ZONE
*/
BOOL StoreUserDisulphide(char *resspec1, char *resspec2)
{
   ZONE *p;
   
   if(gUserSS == NULL)
   {
      INIT(gUserSS, ZONE);
      p = gUserSS;
   }
   else
   {
      p = gUserSS;
      LAST(p);
      ALLOCNEXT(p, ZONE);
   }
   
   if(p==NULL)
   {
      StoreError("StoreUserDisulphide()",
                 "No memory to store data from DISULPHIDE command");
      return(FALSE);
   }
   
   strncpy(p->res1, resspec1, 16);
   strncpy(p->res2, resspec2, 16);

   p->res1[15] = '\0';
   p->res2[15] = '\0';

   return(TRUE);
}

#endif

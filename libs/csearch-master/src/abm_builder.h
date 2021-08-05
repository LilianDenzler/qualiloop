#ifndef __ABM_BUILDER_H__
#define __ABM_BUILDER_H__
/*****************************************************************************
 *      Name: abm_builder.h                                                  *
 *  Function: To define any parameters to be used in the AbM antibody builder*
 *            module.                                                        *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 15/06/93                                                       *
 *---------------------------------------------------------------------------*
 *    Inputs:                                                                *
 *   Outputs:                                                                *
 *   Returns:                                                                *
 * Externals:                                                                *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/

#define MAX_LINE_LEN 240

/* TRUE/FALSE flags if not already set */

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* Batch or interactive run mode definitions */

#define INTERACTIVE  1
#define BATCH        0

/* Some useful floats */

#define ONE          (float)1.0
#define TWO          (float)2.0
#define ZERO         (float)0.0
#define TEN          (float)10.0

/* Some fixed file names */

#define AbM_PROGRESS "PROGRESS.AbM"      /* AbM processing progress file */
#define AbM_SUCCESS  "SUCCESS.AbM"       /* AbM success file */
#define AbM_ERROR    "ERRORS.AbM"        /* The AbM module error output file */
#define FRAME_FILE   "coords.rdb"        /* Framework file from framebldr */
#define CHOTHIA_OUT  "CanonicalReport"   /* CHOTHIA output report file */
#define CONGEN_OUT   "reference.pdb"     /* CONGEN output cooridnates file */
#define EUR_RENAME   "EurRename"         /* Eureka file renaming file name */

/* Step type definitions */

#define NULLSTEP     0   /* NULL step state */
#define FRAMEBUILD   1   /* Frame builder step type identifier */
#define PDBCHECK     2   /* PDB file input step type identifier */
#define SIDECHAIN    3   /* Side chain addition step type identifier */
#define DATABASE     4   /* Database search step type identifier */
#define CONGEN       5   /* CONGEN run step type identifier */
#define EUREKA       6   /* EUREKA run step type identifier */
#define SELECT       7   /* Conformational selection step type identifier */
#define FINISH       8   /* Cleanup step type identifier */

/* Residue build status flags */

#define RES_NOT_BUILT  0  /* Used to mark a residue as not built */
#define RES_BBONE_DONE 1  /* Used to mark residues as having built backbone */
#define RES_ALL_DONE   2  /* Used to mark residues as fully built */
#define IGNORE_A       3  /* Ignore All atoms in a residue */
#define IGNORE_S       4  /* Ignore sidechain atoms in a residue */
#define SCHAIN_ON      5  /* Do sidechain construction */

/* Some CONGEN default parameters */

#define ABM_EPS     50.0
#define ABM_CUTNB    5.0
#define ABM_CUTHB    4.5
#define ABM_CUTHA   90.0

/* Miscellaneous */

#define COMMENT     '!'  /* A file comment line character */
#define MAX_RESOLN  10   /* The maximum number of SDR resolution values */
#define CLOSE_LEN    3   /* The length of a chain closure region */
#define BUILD_LEN    5   /* Minimum length of a residue rebuild region */
#define DBASE_LEN    4   /* Minimum length of a database search region */
#define DELETION    '-'  /* Name of a deletion residue (1 letter code) */
#define MAX_AMINO   20   /* Maximum number of amino acids */
#define MAX_RES_ATM 20   /* Maximum number of atoms per amino acid */

/* Type definitions */

typedef struct STEP STEP;               /* Processing step data structure */
typedef struct ResRange ResRange;       /* Residue ranges structure */
typedef struct FileList FileList;       /* File list data structure */
typedef struct PDBinput PDBinput;       /* PDB file input structure */
typedef struct FinalData FinalData;     /* Final processing instruction data */
typedef struct Constraints Constraints; /* Database search constraints data */
typedef struct LoopData LoopData;       /* Chain fragment rebuild data */
typedef struct DataBase DataBase;       /* Data base search data structure */
typedef struct CongenData CongenData;   /* CONGEN data structure */
typedef struct EurekaData EurekaData;   /* Eureka data structure */
typedef struct Select Select;           /* Final conformation selection data */

/* Structure definitions */

struct STEP
{
   STEP  *next;                     /* Pointer to the next step */
   STEP  *last;                     /* Pointer to the previous step */
   short  StepType;                 /* The type of step this is */
   void  *DataPtr;                  /* Pointer to the step data */
   short  Qactive;                  /* Is the step Active flag */
   int    StepNum;                  /* The serial number of this step */
   char   StepTitle[MAX_LINE_LEN];  /* A description of this step */
};

struct ResRange
{
   ResRange *next;    /* Pointer to the next region */
   int       first;   /* The first residue number in the range */
   int       last;    /* The last residue number in the range */
   int       level;   /* Amount of the residue to ignore (all or sidechain) */
};

struct FileList
{
   FileList *next;                   /* Pointer to the next in the list */
   char      FileName[MAX_LINE_LEN]; /* The file name */
};

struct PDBinput
{
   char PDBfile[MAX_LINE_LEN+1];  /* The input PDB file name */
};

struct FinalData
{
   int  QCleanUp;                 /* Should the work files be removed */
   int  QCleanAll;                /* Should special work files be removed */
   int  QRename;                  /* Should the PDB model file be renamed */
   int  QStdPDB;                  /* Should the termini residues be deleted */
   char PDBfile[MAX_LINE_LEN+1];  /* The final model PDB file name */
};

struct Constraints
{
   Constraints *next;        /* Pointer to the next Constraints structure */
   char         Type[3];     /* The type of constraint (DM or DP) */
   int          Range;       /* The C-alpha range */
   float        MinDis;      /* The minimum distance to the Range C-alpha */
   float        MaxDis;      /* The maximum distance to the Range C-alpha */
};

struct LoopData
{
   DataBase   *DBasePtr;    /* Pointer to database data structure */
   int         QDBase;      /* Indicates if database search required */
   CongenData *CGENPtr;     /* Pointer to CONGEN data structure */
   int         QCGEN;       /* Indicates if CONGEN search required */
   EurekaData *EurPtr;      /* Pointer to EUREKA data structure */
   int         QEureka;     /* Indicates if EUREKA run required */
   Select     *SelPtr;      /* Pointer to final loop selection data */
   int         QSelect;     /* Indicates if final selection is required */
};

struct DataBase
{
   int          QResSpan;               /* Has ResSpan[2] been specified */
   int          ResSpan[2];             /* The first and last loop residues */
   int          QClust1;                /* Has Clust1 been specified */
   float        Clust1Res;              /* Initial cluster resolution */
   int          QClust2;                /* Has Clust2 been specified */
   float        Clust2Res;              /* Secondary cluster resolution */
   int          QNkeep;                 /* Has Nkeep been specified */
   int          Nkeep;                  /* Number of loops to keep */
   int          QCadData;               /* Has CadData been specified */
   char         CadData[MAX_LINE_LEN+1];/* The C-alpha database file name */
   int          QConst;                 /* Have any constraints been given */
   Constraints *ConstPtr;               /* Pointer to the constraints list */
   LoopData    *LoopPtr;                /* Pointer to loop processing info */
   char         FragID[3];              /* The rebuild fragment ID code */
};

struct CongenData
{
   ResRange *TopRange;      /* Pointer to the first sidechain range */
   ResRange *TopIgnore;     /* Pointer to the first IGNORE region */
   int       Qgemax;        /* Has gemax been specified */
   float     gemax;         /* Glycine backbone VDW energy maximum */
   int       Qaemax;        /* Has aemax been specified */
   float     aemax;         /* Alanine backbone VDW energy maximum */
   int       Qpemax;        /* Has pemax been specified */
   float     pemax;         /* Proline backbone VDW energy maximum */
   int       Qpremax;       /* Has premax been specified */
   float     premax;        /* Proline ring energy maximum */
   int       QBackParm;     /* Has BackParm been specified */
   float     BackParm[2];   /* Backbone grid step and VDW maximum */
   int       QSdChParm;     /* Has SdChParm been specified */
   float     SdChParm[2];   /* Sidechain grid step and VDW maximum */
   int       QDataBase;     /* Are there database loops to consider? */
   int       QRebuild;      /* Is there a rebuild region */
   int       ResSpan[2];    /* The first and last rebuild loop residues */
   int       QClosMax;      /* Has ClosMax been specified */
   float     ClosMax;       /* Maximum chain closure VDW energy */
   int       QClose;        /* Is there a closure region */
   int       ClsSpan[2];    /* The first and last closure residues */
   int       QAddHyd;       /* Flag to say whether to add hydrogens */
   int       QFormat;       /* Flag to say whether to build or format a model*/
   LoopData *LoopPtr;       /* Pointer to loop processing info */
   char      FragID[3];     /* The rebuild fragment ID code */
};

struct EurekaData
{
   int       Qmorse;      /* Use Morse or harmonic potential */
   int       Qcross;      /* Use cross terms flag */
   int       QOOplane;    /* Use out of plane terms flag */
   int       Q14Dist;     /* Use 1-4 terms flag */
   int       QUseCut;     /* Use a distance cutoff flag */
   int       QCutOff;     /* Has CutOff been specified? */
   float     CutOff;      /* The cutoff distance if one is applied */
   int       QSmooth;     /* Has Smooth been specified? */
   float     Smooth;      /* The distance cutoff smoothing range */
   int       QNcycle;     /* Has Ncycle been specified */
   int       Ncycle;      /* The number of minimisation cycles */
   int       QNkeep;      /* Has NKeep been specified? */
   int       Nkeep;       /* The number of low energy loops to keep */
   LoopData *LoopPtr;     /* Pointer to loop processing info */
   char      FragID[3];   /* The rebuild fragment ID code */
};

struct Select
{
   int       Qfilter;            /* Use SDR filter or lowest energy */
   float     Resoln[MAX_RESOLN]; /* Angular resolution values for SDR */
   int       NResoln;            /* The number of angular resolution values */
   LoopData *LoopPtr;            /* Pointer to loop processing info */
   char      FragID[3];          /* The rebuild fragment ID code */
};

#endif

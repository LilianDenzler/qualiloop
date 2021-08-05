#ifndef __ABM_MAIN_H__
#define __ABM_MAIN_H__
/*****************************************************************************
 *      Name: abm_main.h                                                     *
 *  Function: To define any parameters to be used in the AbM antibody modell-*
 *            ing program.                                                   *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 06/09/93                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/

#include "abm_pdbio.h"   /* PDB read/write definitions */

#define MAX_TEXT_LEN   240 /* General maximum text length */
#define MAX_UDB_FILES   40 /* Maximum UDB file total, limited by menu size */
#define MAX_SEQ_LEN    300 /* Residues in a chain */
#define MIN_LOOP_LEN     3 /* The minimum allowable loop length */
#define MAX_LOOP_LEN    20 /* The minimum allowable loop length */
#define MIN_FRAME_LEN    5 /* The minimum allowable framework length */
#define MAX_DIR_LEN   1024 /* Maximum directory name length */

/* TRUE/FALSE flags if not already set */

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* Some useful floats and doubles */

#define ONE          (float)1.0
#define TWO          (float)2.0
#define FIVE         (float)5.0
#define ZERO         (float)0.0
#define DZERO       (double)0.0
#define TEN          (float)10.0
#define FIFTEEN      (float)15.0
#define THIRTY       (float)30.0
#define FORTYFIVE    (float)45.0
#define SIXTY        (float)60.0
#define NINETY       (float)90.0
#define ONE_HUNDRED  (float)100.0
#define ONE_TWENTY   (float)120.0
#define ONE_EIGHTY   (float)180.0

/* Highlight selection video attributes */

#define OFF_VIDEO     0
#define DIM_VIDEO     1
#define BOLD_VIDEO    2
#define INVERSE_VIDEO 3
#define UNDERLINE     4
#define DIM_INVERSE   5
#define BOLD_INVERSE  6

/* Process work file tidy up levels */

#define TIDY_NONE 0
#define TIDY_MOST 1
#define TIDY_ALL  2

/* AbM Builder run mode */

#define INTERACTIVE 0
#define BACKGROUND  1

/* Menu ID definitions */

#define MAIN_MENU      0
#define FILEIO_MENU    1
#define BUILDER_MENU   2
#define ALIGN_MENU     3
#define DATABASE_MENU  4
#define SYSTEM_MENU    5
#define FRAME_MENU     6
#define LOOP_MENU      7
#define REBUILD_MENU   8
#define CLOSURE_MENU   9
#define INTERACT_MENU 10
#define DBSRCH_MENU   11
#define CONS_MENU     12
#define CGEN_MENU     13
#define SCHAIN_MENU   14
#define ENERGY_MENU   15

/* Some general menu definitions */

#define SCREEN_WIDTH       80   /* The menu screen width */
#define DATA_LINES         16   /* Lines in the data area */
#define DATA_AREA_TOP       2   /* Line starting the data area */
#define MSG_LINES           4   /* Lines in the message area */
#define MSG_AREA_TOP       18   /* Line starting the message area */
#define PROMPT_LINE        22   /* The line in which to display prompt */
#define PROMPT_TEXT      "=>"   /* The prompt text itself */
#define PROMPT_LEN          3   /* The length of the prompt text + 1 */
#define DELETE              8   /* Delete key character code */

/* Loop construction definitions */

#define REBUILD_LEN      5    /* Length of CONGEN rebuild range */
#define CLOSURE_LEN      3    /* Length of CONGEN closure range */
#define BOND_MORSE       1    /* Bond stretch MORSE type function */
#define BOND_HARMONIC    2    /* Bond stretch harmonic type function */
#define SELECT_EMIN      1    /* Select loop energy loop */
#define SELECT_SDR       2    /* Use SDR ranking to select loop */
#define BLD_CONGEN       1    /* Loop build using CONGEN only */
#define BLD_DATABASE     2    /* Loop build using database only */
#define BLD_COMBINED     3    /* Loop build using combined method */
#define BLD_SIDECHAIN    4    /* Build loop sidechains only */
#define BLD_CANONICAL    5    /* Loop build using canonical method */
#define BLD_NONE         6    /* Loop build using no method */
#define BLD_DEFAULT      7    /* Use the default build method */
#define BLD_INC_NONE     0    /* Include no residue atoms in calculations */
#define BLD_INC_BACK     1    /* Include backbone atoms in calculations */
#define BLD_INC_ALL      2    /* Include all atoms in calculations */
#define BLD_VOID_INC    -1    /* Indeterminable region in build to include*/
#define BLD_VOID_EXC    -2    /* Indeterminable region in build to exclude*/
#define BLD_SIDE_FALSE  10    /* Sidechain construction will not be done */
#define BLD_SIDE_TRUE   11    /* Sidechain construction will be done */
#define BLD_SIDE_VOID  -10    /* Sidechain construction is indeterminable */

#define REBUILD_LEFT_LEFT   10 /* Code for lowering left rebuild limit */
#define REBUILD_LEFT_RIGHT  20 /* Code for raising left rebuild limit */
#define REBUILD_RIGHT_LEFT  30 /* Code for lowering right rebuild limit */
#define REBUILD_RIGHT_RIGHT 40 /* Code for raising right rebuild limit */
#define CLOSURE_LEFT        10 /* Code for moving closure region left */
#define CLOSURE_RIGHT       20 /* Code for moving closure region right */

#define RESIDUE_INTERACT    10 /* Code for residue menu INTERACT mode */
#define RESIDUE_SIDE        20 /* Code for residue menu SIDECHAIN mode */

/* Sequence reading problem definitions */

#define CHAINS_OKAY         30 /* No problem reading chain sequences*/
#define CHAIN_NO_HEADER     31 /* No chain header founf in PDB file */
#define CHAIN_TOO_LONG      32 /* A chain was too long and was truncated */
#define CHAIN_EXCESS        33 /* Too many chains read in */
#define CHAIN_MISSING       34 /* Not enough chains read in */

/* Miscellaneous */

#define NULL_INPUT              100  /* Code for blank user input */
#define REAL_INPUT              101  /* Code for float user input */
#define INT_INPUT               102  /* Code for int user input */
#define CHAR_INPUT              103  /* Code for char user input */
#define VOID_CONSTRAINT (float)-1.0  /* An invalid constraint distance */

/* Type definitions */

typedef struct MenuData  MenuData;        /* Structure for holding menu data */
typedef struct GroupData GroupData;       /* Structure for highlight groups */
typedef struct SysData   SysData;         /* Structure for system data */
typedef struct UDBdata   UDBdata;         /* Structure for UDB file data */
typedef struct SeqData   SeqData;         /* Structure for sequence data */
typedef struct LoopData  LoopData;        /* Structure for loop data */
typedef struct IntrData  IntrData;        /* Interaction menu data structure */
typedef char   ResHead[4][MAX_SEQ_LEN];   /* Interaction menu chain header */
typedef char   ScrnTxt[SCREEN_WIDTH+1];   /* Screen width text string */
typedef int    IncRes[MAX_SEQ_LEN];       /* Array of residues to include */
typedef int    SideRes[MAX_SEQ_LEN];      /* Array of sidechains to build */

/* Structure definitions */

struct MenuData
{
   char       DataLine[DATA_LINES][SCREEN_WIDTH];  /* Data area text */
   char       MsgLine[MSG_LINES][SCREEN_WIDTH];    /* Message area text */
   char       TitleLine[SCREEN_WIDTH];     /* Room for the title line text */
   int        Ngroups;                     /* The number of highlight groups */
   GroupData *Groups;                      /* Pointer to groups array */
};

struct GroupData
{
   int       Nareas;     /* The number of highlight areas in this group */
   int      *Top;        /* Pointer to the array of top most line nos. */
   int      *Bottom;     /* Pointer to the array of bottom most line nos. */
   int      *Left;       /* Pointer to the array of left most column nos. */
   int      *Right;      /* Pointer to the array of right most column nos. */
   ScrnTxt  *KeyWords;   /* Pointer to keyword array associated with areas */
};

struct SysData
{
   /* General menu wise information */

   char CurrentDir[MAX_DIR_LEN]; /* The current working directory */
   int  QLogFile;                /* Whether a log file is specified */
   char LogFile[MAX_TEXT_LEN];   /* The log file name */
   int  RunMode;                 /* The builder process runmode */
   int  QDBname;                 /* Whether user has given a CAD file name */
   char DBname[MAX_TEXT_LEN];    /* CAD file name given by the user */
   int  QPDBdir;                 /* Whether a PDB directory is given */
   char PDBdir[MAX_TEXT_LEN];    /* The PDB directory name */
   char PDBpref[MAX_TEXT_LEN];   /* The PDB prefix name */
   char PDBext[MAX_TEXT_LEN];    /* The PDB file extension */
   int  QExclude;                /* Whether a PDB exclusion list is given */
   char Exclude[MAX_TEXT_LEN];   /* The PDB file exclusion list */
   int  TidyLevel;               /* The process file cleanup level */
   int  HighLite;                /* The highlite/select video attribute */
   int  QReModel;                /* Is the process a remodelling one ? */
   char PDBmodel[MAX_TEXT_LEN];  /* The PDB file to be remodelled */
   int  QModelName;              /* Whether a model name is given */
   char ModelName[MAX_TEXT_LEN]; /* The chosen model name */
   int  MenuType;                /* Menu ID code */
   char PirTitle[MAX_TEXT_LEN];  /* The title for output sequence files */
   int  QSetBuild;               /* Whether the model is ready to build */
   int  QSingleWin;              /* TRUE if MUST NOT use multiple windows */

   /* The specific antibody definition wise information */

   int       Nchains;  /* Copy of the number of configuration chains */
   int       CurChain; /* Current active chain identifier */
   int       Nloops;   /* Copy of the number of configuration loops */
   UDBdata  *UDBptr;   /* Pointer to array of UDB structures */
   SeqData  *SeqPtr;   /* Pointer to array of sequence information */
   LoopData *LoopPtr;  /* Pointer to array of loops data */
   LoopData *CurLoop;  /* Pointer to current active loop */
};

struct UDBdata
{
   int   ChainNum;                          /* The chain number */
   int   NUDBfiles;                         /* The number of UDB files */
   char  Seqs[MAX_UDB_FILES][MAX_SEQ_LEN];  /* Array of sequences */
   float XCalf[MAX_UDB_FILES][MAX_SEQ_LEN]; /* C-alpha X coord array */
   float YCalf[MAX_UDB_FILES][MAX_SEQ_LEN]; /* C-alpha Y coord array */
   float ZCalf[MAX_UDB_FILES][MAX_SEQ_LEN]; /* C-alpha Z coord array */
   char  Names[MAX_UDB_FILES][MAX_TEXT_LEN];/* Array of file names */
   int   QInClass[MAX_UDB_FILES];           /* Is a UDB loop in a length class*/
   float Scores[MAX_UDB_FILES];             /* Array of current scores */
   int   QGetBest;                          /* Choose best score or not */
   int   Chosen;                            /* The current chosen file */
};

struct SeqData
{ 
   int   ChainNum;              /* The chain number for this sequence */
   int   MaxLen;                /* Expected maximum chain length */
   int   Qalign;                /* Whether or not the chain is aligned */
   char  UserSeq[MAX_SEQ_LEN];  /* The users input sequence */
   char  FrameSeq[MAX_SEQ_LEN]; /* The CHOOSER version of UserSeq */
   int   Nfrags;                /* Copy of the number of fragments */
   int  *FrgStart;              /* Array of fragment start residues */
   int  *FrgEnd;                /* Array of fragment end residue */
};

struct LoopData
{
   /* General information */

   int      ChainNum;          /* Chain serial number containing it */
   int      FragNum;           /* Fragment number in ChainNumth chain */
   char     Name[MAX_TEXT_LEN];/* The name of this loop */
   int      Length;            /* The loop length */
   int      BuildMethod;       /* The loop construction method */
   int      Priority;          /* The construction priority */
   int      QCanonical;        /* Does it have a canonical class...*/
   int      CanClass;          /* ...if so the class number */
   IncRes  *Includes;          /* Array of inclusion residues (1 per chain)*/
   SideRes *SideChain;         /* Array of sidechain builds (1 per chain)*/
 
   /* Database search information */

   float InitRes;              /* First cluster resolution */
   int   NDBkeep;              /* Secondary cluster target */
   float MaxRes;               /* Maximum secondary cluster resolution */
   float DPmin[MAX_LOOP_LEN];  /* DP minimum constraint distances */
   float DPmax[MAX_LOOP_LEN];  /* DP maximum constraint distances */
   float DMmin[MAX_LOOP_LEN];  /* DM minimum constraint distances */
   float DMmax[MAX_LOOP_LEN];  /* DM maximum constraint distances */

   /* CONGEN parameters */

   float gemax;      /* Glycine backbone energy maximum */
   float aemax;      /* Alanine backbone energy maximum */
   float pemax;      /* Proline backbone energy maximum */
   float premax;     /* Proline ring energy maximum */
   float bgrid;      /* Backbone search grid size */
   float bmvdw;      /* Backbone VDW energy limit */
   float cmvdw;      /* Chain closure VDW limit */
   float sgrid;      /* Sidechain search grid size */
   float smvdw;      /* Sidechain VDW energy limit */
   int   Qrebuild;   /* Is there a rebuild region or just a closure */
   int   Rebuild[2]; /* Start + stop rebuild residues, wrt loop start */
   int   Closure[2]; /* Start + stop closure residues, wrt loop start  */

   /* Eureka parameters */

   int   Qcross;     /* Use cross terms or not */
   int   Qooplane;   /* Use out of plane terms or not */
   int   Q14dis;     /* Use 1-4 non-bond interactions or not */
   int   Qcutoff;    /* Use a distance cutoff or not */
   float DistCut;    /* The distance cutoff */
   float Smooth;     /* The cutoff smoothing range */
   int   BndModel;   /* The bond stretch model to use */
   int   Ncycle;     /* The number of minimisation cycles to use */
   int   NEkeep;     /* Number of low energy loop models to keep */
   
   /* The final loop selection parameters */

   int   QFMethod;    /* The loop selection method to use */
   float SDRinterval; /* The interval between resolution attempts */
   int   Nintervals;  /* The number of attempts to make */
};

struct IntrData
{
   /* Textural display area variables */

   char       Title[SCREEN_WIDTH+1];                  /* Title line */
   char       DataLines[DATA_LINES][SCREEN_WIDTH+1];  /* Data area text */
   char       MsgLines[MSG_LINES][SCREEN_WIDTH+1];    /* Message area text */

   /* Other data */

   int        Nchains;     /* Number of chains */
   ResHead   *ChainHdr;    /* Chain scale and limit bars */
   int       *ScreenStart; /* Screen left margin w.r.t. chain sequence */
   int        EditMode;    /* The screen edit mode */
   int        EditStart;   /* The edit region first residue */
   int        EditStat;    /* The edit status */
   int        CurChain;    /* The current active chain number */
   int        CurPos;      /* The cursor position on the screen */
   int        QReDraw;     /* Whether to redraw the screen or not */
};

#endif

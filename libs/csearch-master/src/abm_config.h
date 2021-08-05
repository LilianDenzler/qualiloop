#ifndef __ABM_CONFIG_H__
#define __ABM_CONFIG_H__
/*****************************************************************************
 *      Name: abm_config.h                                                   *
 *  Function: Variable definitions for the AbM configuration file reader     *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 03/11/93                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/

#define MAX_STR_LEN   255                /* Maximum string length */
#define MAX_CHAIN_LEN 300                /* Maximum chain length */
#define OML_ABM_ENV   "OML_ABM"          /* Pointer to ABM installation */
#define AbM_CONFIG    "abm_config.dat"   /* The AbM configuration file */

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

typedef struct FragDef                  /* Chain fragment descriptors */
{
   int        FragNum;                  /* The serial number of this fragment */
   int        ChainNum;                 /* The chain this fragment is in */
   int        QIsLoop;                  /* Loop or framework */
   int        CanFlag;                  /* Does it have canonical classes */
   int        UseCan;                   /* Whether or not to use it */
   int        Class;                    /* Option store of class number */
   int        start;                    /* The first residue number */
   int        finish;                   /* The last residue number */
   char       Name[MAX_STR_LEN+1];      /* The fragment name */
} FRAGCONF;

typedef struct ChainDef                 /* Chain descriptors */
{
   char      ChainName[MAX_STR_LEN+1];  /* Chain name */
   char      ChainID;                   /* Single character ID */
   int       ChainNum;                  /* The serial number of this chain */
   int       start;                     /* Starting residue, normalised  */
   int       finish;                    /* Last residue, normalised  */
   int       Nfrags;                    /* The number of fragments */
   FRAGCONF *FragData;                  /* Pointer to array of fragment data */
   char      Sequence[MAX_CHAIN_LEN+1]; /* Optional chain sequence store */
   int       QBuilt[MAX_CHAIN_LEN];     /* Indicates residue build completion*/
   char      UDBname[MAX_STR_LEN+1];    /* Optional UDB structure name */
   int       Qdefined;                  /* Optional chain defined flag */
} CHAINCONF;

typedef struct ConfDef                   /* Configuration information */
{
   int        Nchains;                   /* The number of chains in the model */
   char       title[MAX_STR_LEN+1];      /* The title of this configuration */
   CHAINCONF *ChainData;                 /* Pointer to array of chain data */
   int        placeflag;                 /* Placement flag used in framebldr */
   char       AminoFile[MAX_STR_LEN+1];  /* Database search amino acid file */
   char       ToplgyFile[MAX_STR_LEN+1]; /* Structure topology file */
   char       ParamFile[MAX_STR_LEN+1];  /* Energy parameters file */
   char       SCtopFile[MAX_STR_LEN+1];  /* Side chain topology file */
   char       AlaMapFile[MAX_STR_LEN+1]; /* Alanine energy map file */
   char       GlyMapFile[MAX_STR_LEN+1]; /* Glycine energy map file */
   char       ProMapFile[MAX_STR_LEN+1]; /* Proline energy map file */
   char       ProConFile[MAX_STR_LEN+1]; /* Proline ring construction file */
   char       PGPFile[MAX_STR_LEN+1];    /* Proton generation parameter file */
   char       EurParms[MAX_STR_LEN+1];   /* Eureka parameters file name */
   char       EurIntTypes[MAX_STR_LEN+1];/* Eureka interface atom types */
   char       MutMatrix[MAX_STR_LEN+1];  /* Dayhoff mutation matrix file */
   char       CADfile[MAX_STR_LEN+1];    /* The C-alpha database file name */
} CONFIG;

#endif

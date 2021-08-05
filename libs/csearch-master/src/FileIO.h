/****************************************************************************** 
 *      Name: FileIO.h                                                        *
 *  Function: Header file - contains definitions of constants for use in the  *
 *            FileIO, dialogs and functions.                                  *
 * Copyright: (C) Oxford Molecular Limited.                                   *
 *----------------------------------------------------------------------------*
 *    Author: K.J. Woods, Tessella Support Services Ltd.                      *
 *      Date: 26/03/90                                                        *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * dd/mm/yy                                                                   *
 * 23/07/90 GCC     Added declarations for directory IO operations            *
 * 04/02/91 GCC     Modified scan minimum time to 0.                          *
 * 26/03/91 GCC     Added new file formats.                                   *
 * 15/04/91 GCC     Added scan molecule number.                               *
 * 16/12/92 KJW     Prevent the DirLink and DirSorted structures from being   *
 *                  included more than once (since they are also defined in   *
 *                  DirInc.h for some unfathomable reason).                   * 
 ******************************************************************************/

#ifndef __FILEIO__
#define __FILEIO__

#define STRING_CSSR       "CSSR"
#define STRING_PDB        "PDB"
#define STRING_LIST       "LIST"
#define STRING_MOL        "MOL"
#define STRING_SD         "SD"
#define STRING_COSMIC     "COSMIC"
#define STRING_EUREKA     "EUREKA"

#define OLD               0            /* flag new molecule */
#define NEW               1            /* flag old molecule */

/* Define the maximum length of a file name string */
#define LENVOLNAME 27
#define MAXFILLEN 255
#define LENFILENAME MAXFILLEN-1
#define LENSAVENAME 14   /* maximum length of save file name */

/* maximum path name length */
#ifndef PATH_MAX
#define PATH_MAX 1023
#endif

/* maximum file name length */
#ifndef NAME_MAX
#ifdef __HPUX_7__
#define NAME_MAX 14
#else
#define NAME_MAX 255
#endif
#endif

/* declare structure for an element of a linked list used to hold
   file names from a directory listing */
#if !defined(__Dir_Link__)
#define __Dir_Link__
typedef struct
   {
   char *BackLink;
   char *ForLink;
   char FileName[MAXFILLEN];
   } DirLink, *DirLinkPtr;
#endif

/* declare the size of an element of this list */
#if !defined(LINKSIZE)
#define LINKSIZE sizeof(DirLink)
#endif

/* declare an array element for a list of sorted file names, and their 
   selection status, selected or not selected */

#if !defined(__Dir_Sorted__)
#define __Dir_Sorted__
typedef struct
   {
   char FileName[MAXFILLEN];
   int SelectedFlag;
   } DirSorted, *DirSortedPtr;
#endif

/* declare constants to state of a fileis selected or not */
enum { DIR_NOT_SELECTED, DIR_SELECTED } ;
enum { REMOVE_MOLECULE_CLEAR,
       REMOVE_MOLECULE_LIST,
       REMOVE_MOLECULE_CANCEL,
       REMOVE_MOLECULE_ENTER } ;

/* define constants for scan utility setup */
#define SCAN_SETUP_MAX_OPTIONS 4

enum {	SCAN_SETUP_ADD_FILES,
	SCAN_SETUP_LIST_FILES,
	SCAN_SETUP_CLEAR_LIST,
	SCAN_SETUP_FINISH
     };

#define SCAN_MAX_OPTIONS 12
enum {	SCAN_FIRST,
	SCAN_PREVIOUS,
	SCAN_NEXT,
	SCAN_LAST,
	SCAN_STRUCTURE,
	SCAN_SKIP,
	SCAN_LOOP,
	SCAN_CYCLE,
	SCAN_STOP,
	SCAN_SET_TIME,
	SCAN_KEEP,
	SCAN_FINISH
     };

/* declare structure for an element of a linked list used to hold
   file names from a directory listing */

typedef struct _ScanLink
   {
   struct _ScanLink *BackLink;
   struct _ScanLink *ForLink;
   char FileName[MAXFILLEN];
   short FileType;
   short Skip;
   short Nstruc;               /* Number of coorindate sets in the file. */
} ScanLink, *ScanLinkPtr;

enum { SCAN_SKIP_STRUCTURE, SCAN_NOSKIP_STRUCTURE };
#define SCAN_NO_STRUCTURE -1

/* the scan application status record */
typedef struct
   {
      ScanLinkPtr   ScanLinkPointer;  /* pointer to the list of structures */
      short         CurrentStructure; /* the 'current' structure */
      short         CurrentCoords;    /* current coords in CurrentStructure */
      short         CurrentMolecule;  /* the 'current' scan molecule number */
      short         CycleDirection;   /* the direction to cycle */
      short         State;            /* state (active or inactive ) */
   } ScanRecord;
 
enum {SCAN_CYCLE_FORWARD, SCAN_CYCLE_BACKWARD};
enum {SCAN_INACTIVE, SCAN_ACTIVE};

#define MAX_STRUCTURE_NAME_LENGTH 15

/* define scan loop/cycle default times */
#define SCAN_DEFAULT_TIME 0
#define SCAN_MIN_TIME 0
#define SCAN_MAX_TIME 3600

#endif

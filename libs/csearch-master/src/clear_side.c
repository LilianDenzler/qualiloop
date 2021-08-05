/*****************************************************************************
 *      Name: clear_side.c                                                   *
 *  Function: To give all the side chain atoms in a certain residue range    *
 *            the coordinates +9999.0 so that the interactions with these    *
 *            atoms will be ignored. Glycine and Proline residues are to be  *
 *            unaffected by this action.                                     *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 08/10/93                                                       *
 *---------------------------------------------------------------------------*
 *    Inputs: res1  - the first residue number in the required range         *
 *            res2  - the last  residue number in the required range         *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CongenProto.h"
 
#define IGNORE_ATOM (float)9999.0

void clear_side(char *res1, char *res2)
{
   int atom_start; /* First atom number in first residue */
   int atom_stop;  /* Last atom number in Last residue */
   int iatom;      /* Used to cycle between atom_start and atom_stop */
   int ires;       /* Used to cycle over residues */
   int FstRes;     /* res1 in integer form */
   int LstRes;     /* res2 in integer form */
 
   if (sscanf(res1, "%d", &FstRes) != 1)
   {
      fprintf(out,"ERROR: SCLEAR argument %s is not an integer!\n",res1);
      die();
   }
   else if (sscanf(res2, "%d", &LstRes) != 1)
   {
      fprintf(out,"ERROR: SCLEAR argument %s is not an integer!\n",res2);
      die();
   }
   else if ((FstRes < 1) || (LstRes > values.nres))
   {
      fprintf(out,"ERROR: SCLEAR residues out of range\n");
      die();
   }
   else
   {
      for (ires = FstRes; ires <= LstRes; ires++)
      {
         if (strncmp(pstruct.resnme[ires],"PRO",3))
         {
            atom_start = (int)pstruct.lstatm[ires]+1;
            atom_stop  = (int)pstruct.lstatm[ires+1];
            for(iatom = atom_start; iatom <= atom_stop; iatom++)
            {
               if (strncmp(pstruct.atmnme[iatom-1],"N   ",4) &&
                   strncmp(pstruct.atmnme[iatom-1],"H   ",4) &&
                   strncmp(pstruct.atmnme[iatom-1],"CA  ",4) &&
                   strncmp(pstruct.atmnme[iatom-1],"CB  ",4) &&
                   strncmp(pstruct.atmnme[iatom-1],"C   ",4) &&
                   strncmp(pstruct.atmnme[iatom-1],"O   ",4))
               {
                  coords.xcart[iatom-1] = IGNORE_ATOM;
                  coords.ycart[iatom-1] = IGNORE_ATOM;
                  coords.zcart[iatom-1] = IGNORE_ATOM;
               }
            }
         }
      }
   }
}

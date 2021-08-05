#ifndef __PARAMS_H__
#define __PARAMS_H__

/***********************************************************************
 *      NAME: PARAMS.INC                                               *
 *  FUNCTION: To declare the main system parameters                    *
 * COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
 *---------------------------------------------------------------------*
 *    AUTHOR: Robert Williams                                          *
 *      DATE: 05/10/92                                                 *
 *---------------------------------------------------------------------*
 *    INPUTS:                                                          *
 *   OUTPUTS:                                                          *
 *    LOCALS:                                                          *
 *   GLOBALS:                                                          *
 *     CALLS:                                                          *
 *---------------------------------------------------------------------*
 * MODIFICATION RECORD                                                 *
 * DD/MM/YY   INITS   COMMENTS                                         *
 ***********************************************************************
   Parameter Name                    Description
   --------------                    -----------
      MAXAT      -    Maximum number of atoms
      MAXBND     -    Maximum number of bonds
      MAXANG     -    Maximum number of angles
      MAXTOR     -    Maximum number of proper torsion angles
      MAXIMP     -    Maximum number of improper torsion angles
      MAXHB      -    Maximum number of hydrogen-bonds
      MAXNB      -    Maximum number of non-bond pair exclusions
      MAXRES     -    Maximum number of residues
      MAXSEG     -    Maximum number of segments
      MXDORA     -    Maximum number of hydrogen-bond donors or acceptors
      MAXIC      -    Maximum number of internal coordinates
      MAXBP      -    Maximum number of bond parameters
      MAXAP      -    Maximum number of bond angle parameters
      MAXPTP     -    Maximum number of proper torsion parameters
      MAXITP     -    Maximum number of improper torsion parameters
      MAXHBP     -    Maximum number of hydrogen bond parameters
      MAXNBP     -    Maximum number of non-bonded pair parameters
      MAXATU     -    Maximum number of atom types in use
      MAXATT     -    Maximum number of possible atom types
      MXCBUF     -    The main command line buffer size

C Declarations */

#define maxat   6150
#define maxbnd  6250
#define maxang  9150
#define maxtor  3600
#define maximp  3250
#define maxhb   4100
#define maxnb   16150
#define maxres  1050
#define maxseg  20
#define mxdora  1200
#define maxic   6150
#define maxbp   150
#define maxap   350
#define maxptp  75
#define maxitp  55
#define maxhbp  250
#define maxatu  40

/* MAXNBP is related to MAXATU by the formula (MAXATU**2+MAXATU)/2 */
#define maxnbp  1640

#define maxatt  100
#define mxcbuf  10000

#endif

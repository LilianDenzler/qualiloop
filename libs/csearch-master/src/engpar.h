/***********************************************************************
 *      NAME: ENGPAR                                                   *
 *  FUNCTION: To declare the parameters involved in energy calculations*
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
 
  Variable  Array bounds                Description
  --------  ------------                -----------
  EQBDIS     MAXBP       List of equilibrium bond lengths
  BNDCON     MAXBP       List of bond force constants
  EQANG      MAXAP       List of equilibrium bond angles
  ANGCON     MAXAP       List of bond angle force constants
  TORPHS     MAXPTP      Torsion expression phase shift term
  TORMLT     MAXPTP      Torsion angle multiplicity
  TORCON     MAXPTP      Torsion angle force constant
  EQITAN     MAXITP      Equilibrium improper torsion angle list
  IMPCON     MAXITP      Improper torsion angle force constants
  VDWR12     MAXNBP      Lennard-Jones R^12 non-bond coefficient
  VDWR6      MAXNBP      Lennard-Jones R^6 non-bond coefficient
  HBR12      MAXHBP      R to the power 12 H-bond energy terms
  HBR10      MAXHBP      R to the power 10 H-bond energy terms
  ATMPOL     MAXATT      Atomic polarizability for each atom
  ATNEFF     MAXATT      Number of outer shell electrons (effective)
  VDWRAD     MAXATT      List of atomic van der Waals radii
  ATFLAG     MAXATT      A flag to identify atoms when read from the
                         main topology file
  BNDKEY     MAXBP       An array to key a particular bond
  ANGKEY     MAXAP       An array to key a particular bond angle
  IMPKEY     MAXITP      An array to key a particular improper torsion
  TORKEY     MAXPTP      An array to key a particular proper torsion
  NBKEY      MAXNBP      An array to key a particular non-bond interaction 
  HBKEY      MAXHBP      An array to key a particular hydrogen-bond
  NBCUT        -         The non-bond cutoff distance
  DIELEC       -         The dielectric constant
  NBFLAG       -         The non-bonded calculation options flag

  Declarations  */

#ifndef __ENGPAR_H__
#define __ENGPAR_H__

extern struct {
       float eqbdis[maxbp], bndcon[maxbp],
             eqang[maxap], angcon[maxap],
             torphs[maxptp], tormlt[maxptp], torcon[maxptp],
             eqitan[maxitp], impcon[maxitp],
             vdwr12[maxnbp], vdwr6[maxnbp],
             hbr12[maxhbp], hbr10[maxhbp],
             atmpol[maxatt],
             atneff[maxatt],
             vdwrad[maxatt];

       short atflag[maxatt];

       int   bndkey[maxbp],
             angkey[maxap],
             impkey[maxitp],
             torkey[maxptp],
             nbkey[maxnbp],
             hbkey[maxhbp],
             nbcut, 
             dielec, 
             nbflag;

       } engpar;

#endif

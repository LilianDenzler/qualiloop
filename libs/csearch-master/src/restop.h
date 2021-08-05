/***********************************************************************
 *      NAME: RESTOP                                                   *
 *  FUNCTION: To declare the variable arrays used to store the         *
 *            information read in from the main residue topology file. *
 *            These are then used to derive the potential energy       *
 *            function parameters                                      *
 * COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
 *---------------------------------------------------------------------*
 *    AUTHOR: Robert Williams                                          *
 *      DATE: 06/10/92                                                 *
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
  Variable   Array Bounds                 Description
  --------   ------------                 -----------

  Parameter variables:

  MXATMR         -        The maximum number of atoms per residue
  MXBNDR         -        The maximum number of bonds per residue
  MXANGR         -        The maximum number of bond angles per residue
  MXTORR         -        The maximum number of torsions per residue
  MXIMPR         -        The maximum number of improper torsions per
                          residue
  MXNBER         -        The maximum number of non-bond exclusions
                          per residue
  MXHBAR         -        The maximum number of H-bond acceptors per
                          residue
  MXHBDR         -        The maximum number of H-bond donors per 
                          residue
  MXICR          -        The maximum number of internal coordinates
                          per residue
  MXGRPR         -        The maximum number of groups per residue
  MAXIMR         -        The maximum allowed number of residues
  NINFO          -        Total number of different parameter types
 
  Atom information variables:

  ATMMAS      (MAXATT)    List of atomic masses, indexed with ACODES
  ATMCHG  (MAXIMR,MXATMR) List of atomic charges 
  ACODES      (MAXATT)    Atom type names used with ARMASS
  ACINDX  (MAXIMR,MXATMR) Index of atom type codes for ACODES
  NAMATM  (MAXIMR,MXATMR) List of IUPAC atom names

  Group information variables:

  FSTGRP  (MAXIMR,MXGRPR) List of first atoms in each group
  GRPNAM  (MAXIMR,MXGRPR) List of group names for each residue
  LSTGRP  (MAXIMR,MXGRPR) List of last atoms in each group

  Residue information variables:

  NRESES         -        The number of residues in use
  NAMRES      (MAXIMR)    List of residue names

  Connectivity information:

  BNDAT1  (MAXIMR,MXBNDR) Atom number 1 in each bond
  BNDAT2  (MAXIMR,MXBNDR) Atom number 2 in each bond
  ANGAT1  (MAXIMR,MXANGR) Atom number 1 in each bond angle
  ANGAT2  (MAXIMR,MXANGR) Atom number 2 in each bond angle
  ANGAT3  (MAXIMR,MXANGR) Atom number 3 in each bond angle
  TORAT1  (MAXIMR,MXTORR) Atom number 1 in each proper torsion
  TORAT2  (MAXIMR,MXTORR) Atom number 2 in each proper torsion
  TORAT3  (MAXIMR,MXTORR) Atom number 3 in each proper torsion
  TORAT4  (MAXIMR,MXTORR) Atom number 4 in each proper torsion
  IMPAT1  (MAXIMR,MXIMPR) Atom number 1 in each improper torsion
  IMPAT2  (MAXIMR,MXIMPR) Atom number 2 in each improper torsion
  IMPAT3  (MAXIMR,MXIMPR) Atom number 3 in each improper torsion
  IMPAT4  (MAXIMR,MXIMPR) Atom number 4 in each improper torsion
  EXCLNB  (MAXIMR,MXNBER) List non-bonded exclusions. This is indexed
                          for each atom with NEXCLD
  NEXCLD  (MAXIMR,MXATMR) The number of non-bond exclusions per atom
  ACPTHB  (MAXIMR,MXHBAR) List of hydrogen bond acceptor atoms
  AAN1HB  (MAXIMR,MXHBAR) First antecedent to ACPTHB
  AAN2HB  (MAXIMR,MXHBAR) Second antecedent to ACPTHB
  DONRHB  (MAXIMR,MXHBDR) List of hydrogen bond donor atoms
  DAN1HB  (MAXIMR,MXHBDR) First antecedent to DONRHB
  DAN2HB  (MAXIMR,MXHBDR) Second antecedent to DONRHB
  DHYDHB  (MAXIMR,MXHBDR) Hydrogens attached to DONRHB in H-bond
  NPARAM  (MAXIMR,NINFO)  The total number of some system parameters.
                          atoms, bonds, angles, torsions, improper
                          torsions, non-bond exclusions, H-bond donors,
                          H-bond acceptors, internal coordinates and
                          groups
  NPARUS      (NINFO)     As NPARAM but this gives the number which are
                          actually to be used.
  NPARMX      (NINFO)     Array bounds for NPARAM

  Build information:

  BLDAT1  (MAXIMR,MXICR)  First atom for each internal coordinate in build
  BLDAT2  (MAXIMR,MXICR)  Second atom for each internal coordinate in build
  BLDAT3  (MAXIMR,MXICR)  Third atom for each internal coordinate in build
  BLDAT4  (MAXIMR,MXICR)  Fourth atom for each internal coordinate in build
  B1BLD   (MAXIMR,MXICR)  First bond length for internal coordinates
  B2BLD   (MAXIMR,MXICR)  Next bond length for internal coordinates
  A1BLD   (MAXIMR,MXICR)  First bond angle for internal coordinates
  A2BLD   (MAXIMR,MXICR)  Next bond angle for internal coordinates
  TORBLD  (MAXIMR,MXICR)  Torsion angle value for internal coordinates
*/

#ifndef __RESTOP_H__
#define __RESTOP_H__

#define mxatmr  70
#define mxbndr 100
#define mxangr 150
#define mxtorr  35
#define mximpr  22
#define mxnber  70
#define mxhbar   8
#define mxhbdr   6
#define mxicr   70
#define mxgrpr  10
#define maximr  40
#define ninfo   10

extern struct {

      short bndat1[maximr][mxbndr], bndat2[maximr][mxbndr],
            angat1[maximr][mxangr], angat2[maximr][mxangr],
            angat3[maximr][mxangr], torat1[maximr][mxtorr],
            torat2[maximr][mxtorr], torat3[maximr][mxtorr],
            torat4[maximr][mxtorr], impat1[maximr][mximpr],
            impat2[maximr][mximpr], impat3[maximr][mximpr],
            impat4[maximr][mximpr], exclnb[maximr][mxnber],
            nexcld[maximr][mxatmr], acpthb[maximr][mxhbar],
            aan1hb[maximr][mxhbar], aan2hb[maximr][mxhbar],
            donrhb[maximr][mxhbar], dan1hb[maximr][mxhbar],
            dan2hb[maximr][mxhbar], dhydhb[maximr][mxhbar],
            fstgrp[maximr][mxgrpr], lstgrp[maximr][mxgrpr],
            acindx[maximr][mxatmr], nparam[maximr][ninfo], 
            nparus[ninfo], nparmx[ninfo];

      int   nreses;

      float atmmas[maxatt], atmchg[maximr][mxatmr],
            b1bld[maximr][mxicr], b2bld[maximr][mxicr],
            a1bld[maximr][mxicr], a2bld[maximr][mxicr],
            torbld[maximr][mxicr];
        
      char  bldat1[maximr][mxicr][4], bldat2[maximr][mxicr][4],
            bldat3[maximr][mxicr][4], bldat4[maximr][mxicr][4], 
            namres[maximr][4], acodes[maxatt][4],
            namatm[maximr][mxatmr][4], 
            grpnam[maximr][mxgrpr][4];

       } restop;

#endif

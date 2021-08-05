C***********************************************************************
C*      NAME: RESTOP                                                   *
C*  FUNCTION: To declare the variable arrays used to store the         *
C*            information read in from the main residue topology file. *
C*            These are then used to derive the potential energy       *
C*            function parameters                                      *
C* COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
C*---------------------------------------------------------------------*
C*    AUTHOR: Robert Williams                                          *
C*      DATE: 06/10/92                                                 *
C*---------------------------------------------------------------------*
C*    INPUTS:                                                          *
C*   OUTPUTS:                                                          *
C*    LOCALS:                                                          *
C*   GLOBALS:                                                          *
C*     CALLS:                                                          *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DD/MM/YY   INITS   COMMENTS                                         *
C***********************************************************************
C Variable   Array Bounds                 Description
C --------   ------------                 -----------

C Parameter variables:

C MXATMR         -        The maximum number of atoms per residue
C MXBNDR         -        The maximum number of bonds per residue
C MXANGR         -        The maximum number of bond angles per residue
C MXTORR         -        The maximum number of torsions per residue
C MXIMPR         -        The maximum number of improper torsions per
C                         residue
C MXNBER         -        The maximum number of non-bond exclusions
C                         per residue
C MXHBAR         -        The maximum number of H-bond acceptors per
C                         residue
C MXHBDR         -        The maximum number of H-bond donors per 
C                         residue
C MXICR          -        The maximum number of internal coordinates
C                         per residue
C MXGRPR         -        The maximum number of groups per residue
C MAXIMR         -        The maximum allowed number of residues
C NINFO          -        Total number of different parameter types
 
C Atom information variables:

C ATMMAS      (MAXATT)    List of atomic masses, indexed with ACODES
C ATMCHG  (MXATMR,MAXIMR) List of atomic charges 
C ACODES      (MAXATT)    Atom type names used with ARMASS
C ACINDX  (MXATMR,MAXIMR) Index of atom type codes for ACODES
C NAMATM  (MXATMR,MAXIMR) List of IUPAC atom names

C Group information variables:

C FSTGRP  (MXGRPR,MAXIMR) List of first atoms in each group
C GRPNAM  (MXGRPR,MAXIMR) List of group names for each residue
C LSTGRP  (MXGRPR,MAXIMR) List of last atoms in each group

C Residue information variables:

C NRESES         -        The number of residues in use
C NAMRES      (MAXIMR)    List of residue names

C Connectivity information:

C BNDAT1  (MXBNDR,MAXIMR) Atom number 1 in each bond
C BNDAT2  (MXBNDR,MAXIMR) Atom number 2 in each bond
C ANGAT1  (MXANGR,MAXIMR) Atom number 1 in each bond angle
C ANGAT2  (MXANGR,MAXIMR) Atom number 2 in each bond angle
C ANGAT3  (MXANGR,MAXIMR) Atom number 3 in each bond angle
C TORAT1  (MXTORR,MAXIMR) Atom number 1 in each proper torsion
C TORAT2  (MXTORR,MAXIMR) Atom number 2 in each proper torsion
C TORAT3  (MXTORR,MAXIMR) Atom number 3 in each proper torsion
C TORAT4  (MXTORR,MAXIMR) Atom number 4 in each proper torsion
C IMPAT1  (MXIMPR,MAXIMR) Atom number 1 in each improper torsion
C IMPAT2  (MXIMPR,MAXIMR) Atom number 2 in each improper torsion
C IMPAT3  (MXIMPR,MAXIMR) Atom number 3 in each improper torsion
C IMPAT4  (MXIMPR,MAXIMR) Atom number 4 in each improper torsion
C EXCLNB  (MXNBER,MAXIMR) List non-bonded exclusions. This is indexed
C                         for each atom with NEXCLD
C NEXCLD  (MXATMR,MAXIMR) The number of non-bond exclusions per atom
C ACPTHB  (MXHBAR,MAXIMR) List of hydrogen bond acceptor atoms
C AAN1HB  (MXHBAR,MAXIMR) First antecedent to ACPTHB
C AAN2HB  (MXHBAR,MAXIMR) Second antecedent to ACPTHB
C DONRHB  (MXHBDR,MAXIMR) List of hydrogen bond donor atoms
C DAN1HB  (MXHBDR,MAXIMR) First antecedent to DONRHB
C DAN2HB  (MXHBDR,MAXIMR) Second antecedent to DONRHB
C DHYDHB  (MXHBDR,MAXIMR) Hydrogens attached to DONRHB in H-bond
C NPARAM  (NINFO,MAXIMR)  The total number of some system parameters.
C                         atoms, bonds, angles, torsions, improper
C                         torsions, non-bond exclusions, H-bond donors,
C                         H-bond acceptors, internal coordinates and
C                         groups
C NPARUS      (NINFO)     As NPARAM but this gives the number which are
C                         actually to be used.
C NPARMX      (NINFO)     Array bounds for NPARAM

C Build information:

C BLDAT1  (MXICR,MAXIMR)  First atom for each internal coordinate in build
C BLDAT2  (MXICR,MAXIMR)  Second atom for each internal coordinate in build
C BLDAT3  (MXICR,MAXIMR)  Third atom for each internal coordinate in build
C BLDAT4  (MXICR,MAXIMR)  Fourth atom for each internal coordinate in build
C B1BLD   (MXICR,MAXIMR)  First bond length for internal coordinates
C B2BLD   (MXICR,MAXIMR)  Next bond length for internal coordinates
C A1BLD   (MXICR,MAXIMR)  First bond angle for internal coordinates
C A2BLD   (MXICR,MAXIMR)  Next bond angle for internal coordinates
C TORBLD  (MXICR,MAXIMR)  Torsion angle value for internal coordinates
C     

C Parameter declarations

      INTEGER MXATMR, MXBNDR, MXANGR, MXTORR, MXIMPR, MXNBER, 
     +        MXHBAR, MXHBDR, MXICR , MXGRPR, MAXIMR, NINFO
    
      PARAMETER (MXATMR = 70, MXBNDR =100, MXANGR =150)
      PARAMETER (MXTORR = 35, MXIMPR = 22, MXNBER =140)
      PARAMETER (MXHBAR =  8, MXHBDR =  6, MXICR  = 70)
      PARAMETER (MXGRPR = 10, MAXIMR = 40, NINFO  = 10)

C Other declarations

      INTEGER*2 BNDAT1(MXBNDR,MAXIMR), BNDAT2(MXBNDR,MAXIMR),
     +          ANGAT1(MXANGR,MAXIMR), ANGAT2(MXANGR,MAXIMR),
     +          ANGAT3(MXANGR,MAXIMR), TORAT1(MXTORR,MAXIMR),
     +          TORAT2(MXTORR,MAXIMR), TORAT3(MXTORR,MAXIMR),
     +          TORAT4(MXTORR,MAXIMR), IMPAT1(MXIMPR,MAXIMR),
     +          IMPAT2(MXIMPR,MAXIMR), IMPAT3(MXIMPR,MAXIMR),
     +          IMPAT4(MXIMPR,MAXIMR), EXCLNB(MXNBER,MAXIMR),
     +          NEXCLD(MXATMR,MAXIMR), ACPTHB(MXHBAR,MAXIMR),
     +          AAN1HB(MXHBAR,MAXIMR), AAN2HB(MXHBAR,MAXIMR),
     +          DONRHB(MXHBAR,MAXIMR), DAN1HB(MXHBAR,MAXIMR),
     +          DAN2HB(MXHBAR,MAXIMR), DHYDHB(MXHBAR,MAXIMR),
     +          FSTGRP(MXGRPR,MAXIMR), LSTGRP(MXGRPR,MAXIMR),
     +          ACINDX(MXATMR,MAXIMR), NPARAM(NINFO,MAXIMR), 
     +          NPARUS(NINFO,MAXIMR), NPARMX(NINFO)

      INTEGER NRESES

      REAL ATMMAS(MAXATT), ATMCHG(MXATMR), ACODES(MAXATT),
     +     NAMATM(MXATMR,MAXIMR), GRPNAM(MXGRPR,MAXIMR),
     +     NAMRES(MAXIMR), BLDAT1(MXICR,MAXIMR),
     +     BLDAT2(MXICR,MAXIMR), BLDAT3(MXICR,MAXIMR),
     +     BLDAT4(MXICR,MAXIMR), B1BLD(MXICR,MAXIMR), 
     +     B2BLD(MXICR,MAXIMR), A1BLD(MXICR,MAXIMR),
     +     A2BLD(MXICR,MAXIMR), TORBLD(MXICR,MAXIMR)

      COMMON /RESTOP/ BNDAT1, BNDAT2, ANGAT1, ANGAT2, ANGAT3,
     +                TORAT1, TORAT2, TORAT3, TORAT4, IMPAT1,
     +                IMPAT2, IMPAT3, IMPAT4, EXCLNB, NEXCLD,
     +                ACPTHB, AAN1HB, AAN2HB, DONRHB, DAN1HB,
     +                DAN2HB, DHYDHB, FSTGRP, LSTGRP, ACINDX,
     +                NPARAM, NPARUS, NPARMX, NRESES, ATMMAS,
     +                ATMCHG, B1BLD,  B2BLD,  A1BLD,  A2BLD, 
     +                TORBLD, BLDAT1, BLDAT2, BLDAT3, BLDAT4,
     +                NAMRES, ACODES, NAMATM, GRPNAM
      SAVE   /RESTOP/

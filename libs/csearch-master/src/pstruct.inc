C***********************************************************************
C*      NAME: PSTRUCT                                                  *
C*  FUNCTION: To declare the variables used for storing information    *
C*            related to the whole protein structure. This includes    *
C*            atomic data and topology data                            *
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
C Variable  Array Bounds                Description
C --------  ------------                -----------
C  ATMNME     MAXAT      List of IUPAC atom names 
C  ATCODE     MAXAT      List of atom type codes
C  ATMASS     MAXAT      List of atomic masses
C  ATCHRG     MAXAT      List of atomic charges
C  QMOVE      MAXAT      Indicates if atoms can move (0=Y, 1=N)
C  ATBND1     MAXBND     List of principal atoms in bonds
C  ATBND2     MAXBND     List of second atoms in bonds
C  ATANG1     MAXANG     List of first atoms in bond angles
C  ATANG2     MAXANG     List of second atoms in bond angles
C  ATANG3     MAXANG     List of third atoms in bond angles
C  ATTOR1     MAXTOR     List of first atoms in proper torsions
C  ATTOR2     MAXTOR     List of second atoms in proper torsions
C  ATTOR3     MAXTOR     List of third atoms in proper torsions
C  ATTOR4     MAXTOR     List of fourth atoms in proper torsions
C  ATIMP1     MAXIMP     List of first atoms in improper torsions 
C  ATIMP2     MAXIMP     List of second atoms in improper torsions
C  ATIMP3     MAXIMP     List of third atoms in improper torsions
C  ATIMP4     MAXIMP     List of fourth atoms in improper torsions
C  NBEXCL     MAXNB      List of non-bond exclusions, indexed by LSTEXC 
C  LSTEXC     MAXAT      Gives the last non-bond exclusion in NBEXCL
C                        for an atom
C  HBACPT     MXDORA     List of hydrogen bond acceptors
C  HBAAN1     MXDORA     First antecedent atom to HBACPT
C  HBAAN2     MXDORA     Second antecedent atom to HBACPT
C  HBDONR     MXDORA     List of heavy atom donors
C  HBDHYD     MXDORA     Hydrogen bond donor atom attached hydrogen
C  HBDAN1     MXDORA     First antecedant atom to HBDONR
C  HBDAN2     MXDORA     Second antecedant atom to HBDONR
C  LSTATM   (MAXRES+1)   LSTATM(I+1) gives the last atom in residue I
C  RESNDX     MAXRES     Index of where residues are in the topology arrays
C  RESNME     MAXRES     List of residue names
C  RESID      MAXRES     List of residue indentifications
C  SEGID      MAXSEG     List of segment indentifications
C  SEGNDX  (10,MAXSEG+1) This array gives information regarding the
C                        separation of the structure into segments.
C                        The first array argument indicates:
C                              1 Total number of residues
C                              2 Total number of atoms
C                              3 Total number of bonds
C                              4 Total number of bond angles
C                              5 Total number of proper torsions
C                              6 Total number of improper torsions
C                              7 Total non-bond exclusions
C                              8 Total H-bond donor atoms
C                              9 Total H-bond acceptor atoms
C                             10 The segment type flag 
C                        The second argument indicates the last value
C                        for the quantity defined by the first argument
C                        ie SEGNDX(3,I+1) gives the last bond number 
C                        for segment I.

      INTEGER*2 ATBND1(MAXBND), ATBND2(MAXBND), ATANG1(MAXANG), 
     +          ATANG2(MAXANG), ATANG3(MAXANG), ATTOR1(MAXTOR),
     +          ATTOR2(MAXTOR), ATTOR3(MAXTOR), ATTOR4(MAXTOR), 
     +          ATIMP1(MAXIMP), ATIMP2(MAXIMP), ATIMP3(MAXIMP),
     +          ATIMP4(MAXIMP), HBACPT(MXDORA), HBAAN1(MXDORA),
     +          HBAAN2(MXDORA), HBDONR(MXDORA), HBDHYD(MXDORA), 
     +          HBDAN1(MXDORA), HBDAN2(MXDORA), LSTATM(MAXRES+1),
     +          RESNDX(MAXRES), NBEXCL(MAXNB), QMOVE(MAXAT),
     +          ATCODE(MAXAT)

      INTEGER LSTEXC(MAXAT), SEGNDX(10,MAXSEG+1)

      CHARACTER*4 SEGID(MAXSEG), RESID(MAXRES), RESNME(MAXRES), 
     +            ATMNME(MAXAT)
      REAL ATCHRG(MAXAT), ATMASS(MAXAT)

      COMMON /PSTRUCT/SEGID, RESID, RESNME, LSTEXC, SEGNDX,
     +                ATMNME, ATCHRG, ATMASS, 
     +                ATBND1, ATBND2, ATANG1, 
     +                ATANG2, ATANG3, ATTOR1,
     +                ATTOR2, ATTOR3, ATTOR4, ATIMP1, ATIMP2, ATIMP3,
     +                ATIMP4, HBACPT, HBAAN1, HBAAN2, HBDONR, HBDHYD,
     +                HBDAN1, HBDAN2, LSTATM, RESNDX, NBEXCL, QMOVE,
     +                ATCODE
      SAVE   /PSTRUCT/

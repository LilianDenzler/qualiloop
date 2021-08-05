C
C     *************************************************************
C     *                                                           *
C     * File: aeincl.inc                                          *
C     *                                                           *
C     * John Woods, Oxford Molecular Ltd, 1991                    *
C     * (c) Oxford Molecular Ltd 1991                             *
C     *                                                           *
C     * Origin: JAW 02/01/92                                      *
C     *                                                           *
C     * This is the Fortran include file that declares the common *
C     * blocks used by the ABMEURINT conversion utility program   *
C     * that converts PDB files into a format suitable for input  *
C     * to EUREKA.                                                *
C     *                                                           *
C     *************************************************************
C
C     Function: declare the common block AEINCL and the variables stored
C     in it.
C
C
C     Modification history:
C
C     Orginated 02/01/92 JAW
C     11/09/92 JAW Added new variables for CR520
C     02/11/93 RKW Added IGNLEV for CR 2014 to allow sidechains to be
C                  selectively ignored. Also split into CHAR/NON_CHAR
C                  common blocks to comply with ANSI standard.
C
C     *************************************************************
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C     This first block of variables holds data to be written to the
C     .def and .cri files
C
      INTEGER MAXATS,MAXRES,MAXBND
      PARAMETER (MAXATS=10000,MAXRES=300,MAXBND=4)
C*     - maximum number of atoms, residues, bonds per atom
      INTEGER NATOMS,NRESID
C*      - actual numbers of bonds and residues read in
      CHARACTER*4 ATMLAB(MAXATS)
C*     - atom labels
      INTEGER IBONDS(MAXATS,MAXBND)
C*     - holds connectivity information
      REAL CHARGE(MAXATS)
C*     - holds partial charge information
      CHARACTER*12 ATTEUR(MAXATS)
C*     - holds EUREKA atom type information
      REAL ATX(MAXATS),ATY(MAXATS),ATZ(MAXATS)
C*     - holds XYZ coordinates
      INTEGER ATRESN(MAXATS)
C*     - holds a record of which residue each atom belongs to,
C        *after* renumbering by AEPROC (so this information is
C        different from that stored in IFIRST and ILAST below,
C        because CTER and NTER residues have been merged with the
C        following or preceeding residue.
      CHARACTER*3 RSLABN(MAXRES)
C*     - this holds the new residue labels, corresponding to the
C        renumbering in ATRESN.
      INTEGER NEWRSN
C*      - holds the modified number of residues after removing CTER
C        and NTER records
      CHARACTER*3 RESNAM(MAXRES)
C*     - holds residue names
      INTEGER IFIRST(MAXRES),ILAST(MAXRES)
C*     - first and last atom in each residue
C
C     The next block of variables holds data read from the standard
C     database file which contains atom labels, partial charges,
C     connectivity and EUREKA atom types for each atom in each
C     type of residue.
C
      INTEGER EMAXAT,EMAXRS
      PARAMETER (EMAXAT=500,EMAXRS=24)
C*      - maximum number of atoms, residues that can be read from the
C         file
      INTEGER EATOMS, ERESID
C*      - actual number of atoms, residues read from the file
      CHARACTER*3 ERESCO(EMAXRS)
C*      - residue codes
      INTEGER EFIRST(EMAXRS),ELAST(EMAXRS)
C*      - index number of first and last atom in each residue
      CHARACTER*4 EATLAB(EMAXAT)
C*      - atom labels
      INTEGER EATREF(EMAXAT)
C*     - atom reference numbers
      REAL ECHARG(EMAXAT)
C*      - partial charge information
      CHARACTER*12 EEURTY(EMAXAT)
C*      - EUREKA atom types
      INTEGER EIBOND(EMAXAT,MAXBND)
C*      - connectivity info
C
C     The next few variables are used to keep track of the number of
C     molecules we have and where they start and end.
C
      INTEGER MAXMOL
      PARAMETER (MAXMOL=10)
C*      - MAXMOL is the maximum number of molecules we are allowed
      INTEGER NMOLES
C*      - NMOLES is the actual number of molecules present
      INTEGER MOLRS1(MAXMOL)
C*      - MOLRS1 is a list of the last residues in each molecule
C
C
C
      INTEGER MXATRS
      PARAMETER (MXATRS=25)
C*     - the maximum number of atoms in any single residue
C
C     OMUPD JAW 11/09/92 The following are new variables for
C     CR520.
C
      INTEGER RESIGL(MAXRES),NRESIG
      LOGICAL IGNLEV(MAXRES)
C*       - RESIGL is a list of residue numbers that are to be
C          ignored, NRESIG is the number of residues to ignore.
C          The numbers refer to the numbering in the original
C          PDB file - residues are renumbered in the .def
C          file due to the deletions. IGNLEV is TRUE if all atoms
C          are to be excluded else only the sidechains are.
      INTEGER ATCONV(MAXATS,2),NATCNV
C*        - ATCONV is an array holding atom equivalences:
C         Atom ATCONV(i,1) in the input PDB file is equivalent
C         to atom ATCONV(i,2) in the output .def file.  If
C         ATCONV(i,2) is zero then atom i has no equivalent in
C         the .def file (e.g. it is in an ignored region of the
C         input structure).   NATCNV is the total
C         number of atoms in this list of equivalences
C         (this will be larger than NATOMS by the number of
C         ignored atoms). 
C
C     Now put all the necessary variables in common blocks and
C     save them...
C
C Character COMMON block...

      COMMON/AEINCC/ ATMLAB, ATTEUR, RSLABN, RESNAM, ERESCO, EATLAB, 
     +               EEURTY
      SAVE  /AEINCC/

C Integer COMMON block...

      COMMON/AEINCI/ NATOMS, NRESID, IBONDS, ATRESN, NEWRSN, IFIRST, 
     +               ILAST,  EATOMS, ERESID, EFIRST, ELAST,  EATREF,
     +               EIBOND, NMOLES, MOLRS1, RESIGL, NRESIG, ATCONV,
     +               NATCNV, IGNLEV
      SAVE  /AEINCI/

C Real COMMON block...

      COMMON/AEINCR/ CHARGE, ATX, ATY, ATZ, ECHARG
      SAVE  /AEINCR/

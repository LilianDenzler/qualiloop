C***********************************************************************
C*      NAME: abm_makedb.inc                                           *
C*  FUNCTION: To declare the MAKEDB PDB file data                      *
C* COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1993                       *
C*---------------------------------------------------------------------*
C*    AUTHOR: Robert Williams                                          *
C*      DATE: 15/03/93                                                 *
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
C Parameter declarations

C MAXAT  = Maximum number of atoms

      INTEGER MAXAT
      PARAMETER (MAXAT = 15000)

C Character declarations
C ----------------------
C ATMSYM - atomic symbol for each atom
C ATMNAM - the name of each atom
C RESNAM - the name of the residue each atom belongs to
C CHAIN  - the name of the chain each atom belongs to
C INSERT - insert code for each atom

      CHARACTER*2 ATMSYM(MAXAT)
      CHARACTER*4 ATMNAM(MAXAT)
      CHARACTER*4 RESNAM(MAXAT)
      CHARACTER*2 CHAIN(MAXAT)
      CHARACTER*2 INSERT(MAXAT)

C Integer declarations
C --------------------
C NATOMS - the number of atoms in the current PDB file
C RESNUM - the residue number for each atom
C ATMNUM - the serial atom number of each atom

      INTEGER NATOMS
      INTEGER RESNUM(MAXAT)
      INTEGER ATMNUM(MAXAT)

C Real declarations
C -----------------
C XCART  - X cartesian coordinate for each atom
C YCART  - Y cartesian coordinate for each atom
C ZCART  - Z cartesian coordinate for each atom

      REAL XCART(MAXAT)
      REAL YCART(MAXAT)
      REAL ZCART(MAXAT)

C Save variables in common blocks depending on type
C -------------------------------------------------

      COMMON /CHINFO/ ATMSYM, ATMNAM, RESNAM, CHAIN, INSERT
      SAVE   /CHINFO/

      COMMON /IINFO/ NATOMS, RESNUM, ATMNUM
      SAVE   /IINFO/

      COMMON /RINFO/ XCART, YCART, ZCART
      SAVE   /RINFO/

C
C     *************************************************************
C     *                                                           *
C     * File: abmeur.inc                                          *
C     *                                                           *
C     * John Woods, Oxford Molecular Ltd, 1991                    *
C     * (c) Oxford Molecular Ltd 1991                             *
C     *                                                           *
C     * Origin: JAW 19/12/91                                      *
C     *                                                           *
C     * This is the Fortran include file that declares the common *
C     * blocks used by EUREKA to hold flags and other information *
C     * that is used specifically when EUREKA is executed as part *
C     * of the Antibody Modelling package.                        *
C     *                                                           *
C     *************************************************************
C
C     Function: declare the common block ABMEUR and the variables stored
C     in it.
C
C     These variables are copies of C global variables that are declared
C     in the EurekaGlobal.h header file and initialised according to
C     keywords specified in the EUREKA contol deck file.
C
C     Modification history:
C
C     Orginated 20/12/91 JAW
C     10/09/92 JAW Added new variables for CR 520
C
C     *************************************************************
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      LOGICAL QABMMD
C*     - QABMMD specifies whether EUREKA is running in "ABM Mode".
C        When EUREKA is run by ABM this will be .TRUE., otherwise
C        it will be .FALSE.
      INTEGER ABMFLM
      PARAMETER (ABMFLM=250)
C*     - ABMFLM specifies the (maximum) length of a filename string.
C        This should be at least as long as the corresponding C string
C        declared in EurekaGlobal.h
      CHARACTER*(ABMFLM) FABMOE, FABMOS, FABMER, FABMIC, FABMPD,
     +         FABMAT
C*     - the above variables hold various filenames passed to EUREKA
C        by ABM via the control file:
C        FABMOE is the EUREKA output file that will hold the
C               list of calculated energies
C        FABMOS is the root filename for EUREKA output of
C               minimized structures
C        FABMER is the EUREKA error output file for communication with
C               ABM
C        FABMIC is the name of the EUREKA input file holding the list
C               of additional conformations
C        FABMPD is the name of the input PDB file (OMUPD JAW 10/09/92)
C        FABMAT is the name of the input atom conversion file
C                    (OMUPD JAW 10/09/92)
      INTEGER ABMCYC
C*      - ABMCYC is the number of cycles of minimization to be performed
C         on each of the input conformations
      INTEGER ABMKEP
C*      - ABMKEP specifies how many of the lowest-energy structures are
C         to be kept. For example, ABMKEP=5 means the five lowest energy
C         structures found are to be kept.  They will be written to files
C         called root_01, root_02, ..., root_05, where root is the
C         filename held in FABMOS.
      INTEGER MINKEP,MAXKEP
      PARAMETER (MINKEP=5,MAXKEP=50)
C*      - the maximum and minimum legal values for ABMKEP (inclusive)
      LOGICAL QSOLVM
C*      - QSOLVM specifies whether a solvent-modified potential is to be
C         used for energy comparisons - .TRUE. indicates that the
C         solvent-modified potential should be used, .FALSE. indicates
C         that the full potential should be used.
      INTEGER ASTERR
C*      - error status code: flags the status of EUREKA
C         upon completion.  Current values are as follows:
C         0  okay
C         1 Missing or invalid keywords in control deck
C         2 Could not read (or open) conformations file
C        99 Other (fatal) error condition
      INTEGER MXVATM
      PARAMETER (MXVATM=2000)
C*     - MXVATM is the maximum number of atoms that can be
C        specified in the conformations file, i.e. the maximum
C        number of atoms in the conformationally-variable
C        part of the protein
      INTEGER NVATOM
C*      - NVATOM is the actual number of atoms whose coordinates
C         are read from the conformations file
      INTEGER ATLIST(MXVATM)
C*      - ATLIST holds a list of which atom numbers belong
C         to the variable part of the protein
      INTEGER ATCONV(MXTATM,2),NATCNV
C*     - OMUPD JAW 10/09/92 Array holding a list of atom number equivalences.
C       Atom number ATCONV(i,1) in the original PDB file is equivalent to
C       atom number ATCONV(i,2) in the working arrays
C       NATCNV is the number of records in this array.
C
C     Now put all the necessary variables in common blocks and
C     save them...
C
      COMMON/ABMEUR/QABMMD,QSOLVM,
     +              FABMOE,FABMOS,FABMER,FABMIC,
     +              ABMCYC,ABMKEP,ASTERR,
     +              NVATOM,ATLIST,
     +              ATCONV,NATCNV,FABMPD,FABMAT
      SAVE/ABMEUR/

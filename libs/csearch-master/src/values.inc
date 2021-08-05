C***********************************************************************
C*      NAME: VALUES                                                   *
C*  FUNCTION: To declare the most frequently used system values        *
C* COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
C*---------------------------------------------------------------------*
C*    AUTHOR: Robert Williams                                          *
C*      DATE: 05/10/92                                                 *
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
C    Variable              Description
C    --------              -----------
C
C     NATOMS        The number of atoms in the system
C     NRES          The number of residues
C     NSEGS         The number of segments
C     NBONDS        The number of bonds in the system
C     NANGS         The number of bond angles in the system
C     NPTORS        The number of proper torsion angles
C     NITORS        The number of improper torsion angles
C     NHBS          The number of hydrogen bonds
C     NNBS          The number of non-bond pairs
C     NDONAT        The number of hydrogen bond donor atoms
C     NACCAT        The number of hydrogen bond acceptor atoms
C     NBPAR         The number of bond parameters
C     NAPAR         The number of angle parameters
C     NPTPAR        The number of proper torsion parameters
C     NITPAR        The number of improper torsion parameters
C     NHBPAR        The number of hydrogen parameters
C     NATYPS        The number of atom types
C     NBAUTO        Flag to generate non-bonded exclusions
C
C     PI            The value of PI!
C     RAD120        120 degrees expressed as radians
C     DTORAD        Degree TO Radian conversion parameter
C     LARGNUM       Largest acceptable real number
C     LARGINT       Largest acceptable integer
C     ANUM          The coordinates of dummy atoms

C Declarations

      INTEGER NATOMS, NRES, NSEGS, NBONDS, NANGS, NPTORS, NITORS, NHBS,
     +        NNBS, NDONAT, NACCAT, NBPAR, NAPAR, NPTPAR, NITPAR,
     +        NHBPAR, NATYPS, NBAUTO

      COMMON /VALUES/ NATOMS, NRES, NSEGS, NBONDS, NANGS, NPTORS, 
     +                NITORS, NHBS, NNBS, NDONAT, NACCAT, NBPAR,
     +                NAPAR, NPTPAR, NITPAR, NHBPAR, NATYPS, NBAUTO
      SAVE   /VALUES/

      DOUBLE PRECISION PI, RAD120, DTORAD, LARGNUM, ANUM

      INTEGER LARGINT

      PARAMETER (PI = 3.141592653589794)
      PARAMETER (RAD120 = 120.0/180.0*PI, DTORAD = PI/180.0)
      PARAMETER (LARGNUM = 1.7E38, LARGINT = 2147483647, ANUM = 9999.0)

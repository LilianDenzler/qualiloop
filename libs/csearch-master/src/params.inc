C***********************************************************************
C*      NAME: PARAMS.INC                                               *
C*  FUNCTION: To declare the main system parameters                    *
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
C  Parameter Name                    Description
C  --------------                    -----------
C     MAXAT      -    Maximum number of atoms
C     MAXBND     -    Maximum number of bonds
C     MAXANG     -    Maximum number of angles
C     MAXTOR     -    Maximum number of proper torsion angles
C     MAXIMP     -    Maximum number of improper torsion angles
C     MAXHB      -    Maximum number of hydrogen-bonds
C     MAXNB      -    Maximum number of non-bond pair exclusions
C     MAXRES     -    Maximum number of residues
C     MAXSEG     -    Maximum number of segments
C     MXDORA     -    Maximum number of hydrogen-bond donors or acceptors
C     MAXIC      -    Maximum number of internal coordinates
C     MAXBP      -    Maximum number of bond parameters
C     MAXAP      -    Maximum number of bond angle parameters
C     MAXPTP     -    Maximum number of proper torsion parameters
C     MAXITP     -    Maximum number of improper torsion parameters
C     MAXHBP     -    Maximum number of hydrogen bond parameters
C     MAXNBP     -    Maximum number of non-bonded pair parameters
C     MAXATU     -    Maximum number of atom types in use
C     MAXATT     -    Maximum number of possible atom types
C     MXCBUF     -    The main command line buffer size

C Declarations

      INTEGER MAXAT, MAXBND, MAXANG, MAXTOR, MAXIMP, MAXHB, MAXNB, 
     +        MAXRES, MAXSEG, MXDORA, MAXIC, MAXBP, MAXAP, MAXPTP,
     +        MAXITP, MAXHBP, MAXNBP, MAXATU, MAXATT, MXCBUF
   
C Parameterisations
   
      PARAMETER (MAXAT  = 6150)
      PARAMETER (MAXBND = 6250)
      PARAMETER (MAXANG = 9150)
      PARAMETER (MAXTOR = 3600)
      PARAMETER (MAXIMP = 3250)
      PARAMETER (MAXHB  = 4100)
      PARAMETER (MAXNB  = 16150)
      PARAMETER (MAXRES = 1050)
      PARAMETER (MAXSEG = 20)
      PARAMETER (MXDORA = 1200)
      PARAMETER (MAXIC  = 6150)
      PARAMETER (MAXBP  = 150)
      PARAMETER (MAXAP  = 350)
      PARAMETER (MAXPTP = 75)
      PARAMETER (MAXITP = 55)
      PARAMETER (MAXHBP = 250)
      PARAMETER (MAXATU = 40)

C MAXNBP is related to MAXATU by the formula (MAXATU**2+MAXATU)/2
      PARAMETER (MAXNBP = 1640)

      PARAMETER (MAXATT = 100)
      PARAMETER (MXCBUF = 10000)

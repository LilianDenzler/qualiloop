C***********************************************************************
C*      NAME: ENGPAR                                                   *
C*  FUNCTION: To declare the parameters involved in energy calculations*
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
C
C Variable  Array bounds                Description
C --------  ------------                -----------
C EQBDIS     MAXBP       List of equilibrium bond lengths
C BNDCON     MAXBP       List of bond force constants
C EQANG      MAXAP       List of equilibrium bond angles
C ANGCON     MAXAP       List of bond angle force constants
C TORPHS     MAXPTP      Torsion expression phase shift term
C TORMLT     MAXPTP      Torsion angle multiplicity
C TORCON     MAXPTP      Torsion angle force constant
C EQITAN     MAXITP      Equilibrium improper torsion angle list
C IMPCON     MAXITP      Improper torsion angle force constants
C VDWR12     MAXNBP      Lennard-Jones R^12 non-bond coefficient
C VDWR6      MAXNBP      Lennard-Jones R^6 non-bond coefficient
C HBR12      MAXHBP      R to the power 12 H-bond energy terms
C HBR10      MAXHBP      R to the power 10 H-bond energy terms
C ATMPOL     MAXATT      Atomic polarizability for each atom
C ATNEFF     MAXATT      Number of outer shell electrons (effective)
C VDWRAD     MAXATT      List of atomic van der Waals radii
C ATFLAG     MAXATT      A flag to identify atoms when read from the
C                        main topology file
C BNDKEY     MAXBP       An array to key a particular bond
C ANGKEY     MAXAP       An array to key a particular bond angle
C IMPKEY     MAXITP      An array to key a particular improper torsion
C TORKEY     MAXPTP      An array to key a particular proper torsion
C NBKEY      MAXNBP      An array to key a particular non-bond interaction 
C HBKEY      MAXHBP      An array to key a particular hydrogen-bond
C NBCUT        -         The non-bond cutoff distance
C DIELEC       -         The dielectric constant
C NBFLAG       -         The non-bonded calculation options flag
 
C Declarations 

      REAL EQBDIS(MAXBP), BNDCON(MAXBP), EQANG(MAXAP), ANGCON(MAXAP),
     +     TORPHS(MAXPTP), TORMLT(MAXPTP), TORCON(MAXPTP), 
     +     EQITAN(MAXITP), IMPCON(MAXITP), VDWR12(MAXNBP), 
     +     VDWR6(MAXNBP), HBR12(MAXHBP), HBR10(MAXHBP), ATMPOL(MAXATT),
     +     ATNEFF(MAXATT), VDWRAD(MAXATT) 
      
      INTEGER*2 ATFLAG(MAXATT) 

      INTEGER BNDKEY(MAXBP), ANGKEY(MAXAP), IMPKEY(MAXITP),
     +        TORKEY(MAXPTP),NBKEY(MAXNBP), HBKEY(MAXHBP),
     +        NBCUT, DIELEC, NBFLAG

      COMMON /ENGPAR/ EQBDIS, BNDCON, EQANG, ANGCON, TORPHS, TORMLT,
     +                TORCON, EQITAN, IMPCON, VDWR12, VDWR6, HBR12,
     +                HBR10, ATMPOL, ATNEFF, VDWRAD, ATFLAG, BNDKEY,
     +                ANGKEY, IMPKEY, TORKEY, NBKEY, HBKEY, NBCUT,
     +                DIELEC, NBFLAG
      SAVE   /ENGPAR/ 

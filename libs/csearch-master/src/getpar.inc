C***********************************************************************
C*      NAME: GETPAR                                                   *
C*  FUNCTION: To point to parameters for bond etc energy calculations  *
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
C  Variable name   Array bounds             Description
C  -------------   ------------             -----------
C
C      IBNDP        MAXBND      Location of bond parameters 
C      IANGP        MAXANG      Location of angle parameters
C      ITORP        MAXTOR      Location of proper torsion parameters
C      IIMPP        MAXIMP      Location of improper torsion parameters
C      IHBP         MAXHB       Location of hydrogen-bond parameters
C
      INTEGER*2 IBNDP(MAXBND), IANGP(MAXANG), ITORP(MAXTOR), 
     +          IIMPP(MAXIMP), IHBP(MAXHB)

      COMMON /GETPAR/ IBNDP, IANGP, ITORP, IIMPP, IHBP
      SAVE   /GETPAR/

C***********************************************************************
C*      NAME: COORDS.INC                                               *
C*  FUNCTION: To declare the coordinate arrays                         *
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
C  Variable name   array bounds           description
C  -------------   ------------           -----------
C    XCART           (MAXAT)     X cartesian coordinate for each atom
C    YCART           (MAXAT)     Y cartesian coordinate for each atom
C    ZCART           (MAXAT)     Z cartesian coordinate for each atom
C
C    XWORK           (MAXAT)     X work space coordinate for each atom
C    YWORK           (MAXAT)     Y work space coordinate for each atom
C    ZWORK           (MAXAT)     Z work space coordinate for each atom

C Declarations

      REAL XCART(MAXAT), YCART(MAXAT), ZCART(MAXAT)
      REAl XWORK(MAXAT), YWORK(MAXAT), ZWORK(MAXAT)

      COMMON /COORDS/ XCART, YCART, ZCART, XWORK, YWORK, ZWORK
      SAVE   /COORDS/

C***********************************************************************
C*      NAME: HBONDS                                                   *
C*  FUNCTION: To declare hydrogen bonding parameters                   *
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
C  Variable   Array Bounds              Description
C  --------   ------------              -----------
C  
C  HBCUT          -          The hydrogen bond distance limit
C  HBACUT         -          The hydrogen bond angle limit
C  HBDON        MAXHB        The heavy atom donor for each hydrogen bond
C  HBACC        MAXHB        The heavy atom acceptor for each hydrogen
C     
C Declarations

      INTEGER*2 HBDON(MAXHB), HBACC(MAXHB)

      REAL HBCUT, HBACUT

      COMMON /HBONDS/ HBDON, HBACC, HBCUT, HBACUT
      SAVE   /HBONDS/

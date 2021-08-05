C***********************************************************************
C*      NAME: DBG                                                      *
C*  FUNCTION: To declare the DeBuG level parameters for CONGEN         *
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

C Declarations

      INTEGER DBG_CGEN, DBG_CLSCHN, DBG_ALLOC, DBG_ALLSTK, DBG_ALLHP

      COMMON /DBG/ DBG_CGEN, DBG_CLSCHN, DBG_ALLOC, DBG_ALLSTK,
     +             DBG_ALLHP
      SAVE   /DBG/

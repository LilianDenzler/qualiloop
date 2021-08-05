C*****************************************************************************
C*      Name: abmerf                                                         *
C*  Function: to wrap and use the system erf function                        *
C* Copyright: (C) OXFORD MOLECULAR LTD, 1992                                 *
C*---------------------------------------------------------------------------*
C*    Author: Brian Jamieson                                                 *
C*      Date: 16/03/92                                                       *
C*---------------------------------------------------------------------------*
C*    Inputs: x                                                             *
C*   Outputs:                                                                *
C*   Returns: abmerf                                                         *
C* Externals:                                                                *
C*---------------------------------------------------------------------------*
C* MODIFICATION RECORD:                                                      *
C* DD/MM/YY   Initials   Comments                                            *
C*****************************************************************************/
C
      FUNCTION ABMERF(X)
C
      IMPLICIT NONE
C
      REAL ABMERF,X
      REAL XTEMP,XX
C
      XX = X
C
      CALL FWERF(XX,XTEMP)
C
      ABMERF = XTEMP
C
      RETURN
      END

**!!MSGDAT.INS
C
C
C -- Created:-
C      18-Oct-1990      JPH    Define parameters for new message
C                              handling routines.
C
C -- Maximum number of messages.
      INTEGER     MSGMAX
      PARAMETER ( MSGMAX = 100 )
C
C -- The messages.
      CHARACTER*100   MESSAG ( MSGMAX )
      COMMON /MESSGS/ MESSAG
C
      SAVE /MESSGS/
C

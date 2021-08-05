**!!MESSAGES.INS
C
C
C -- Created:-
C      18-Oct-1990      JPH    Define parameters for new message
C                              handling routines.
C
C -- Parameters for 'ZERRM%'.
C
C -- Information message.
      INTEGER MSGTXT
C -- Warning message.
      INTEGER MSGWRN
C -- Error message.
      INTEGER MSGERR
C -- Information message.
      INTEGER MSGINF
C -- Continuation message.
      INTEGER MSGCON
      PARAMETER ( MSGTXT =  0 , MSGWRN =  1 , MSGERR =  2 ,
     +                          MSGINF = -1 , MSGCON = -2 )
C

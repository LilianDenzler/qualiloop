CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MERFI(P,X,IER)
C                **************
C     Calculates the inverse Error Function by doing a binary search
C     of calls to the error function. It finds the top and bottom values
C     which give the same value of P and takes the average of the two.
C     Andrew C.R. Martin     LMB, Oxford.
C     This code may be copied and used by anyone providing this notice
C     is retained
 
      IMPLICIT NONE
 
      INTEGER IER
      REAL SMALL,X,P,BOT,TOP,CUT,ABMERF,ABMERFC,OLD,HI
 
      PARAMETER (SMALL = 0.3E-7)
      LOGICAL START
 
      IF (P .EQ. 0.0) THEN
         X = 0.0
         IER = 0
         RETURN
      ENDIF
 
      IF ((P .GT. 1.0).OR.(P .LT. -1.0)) THEN
         PRINT *, 'Error in MERFI: P out of range'
         PRINT *, 'P = ', P
         IER = 129
         X = 0.99E25
         RETURN
      ENDIF
 
      IF (P.LT.0) THEN
         P = -P
      ENDIF
 
 
C     First search for the maximum value of X which gives this value for P
      BOT = 0.0
      TOP = 5.0
 
      START = .TRUE.
 
  100 CONTINUE
      CUT = BOT + (TOP-BOT)/2
      ABMERFC =ABMERF(CUT)
 
      IF (START) THEN
         START = .FALSE.
      ELSEIF (CUT .EQ. OLD) THEN
         GOTO 200
      ENDIF
 
      OLD = CUT
 
 
      IF (ABMERFC.LE.P) THEN
         IF (ABMERF(CUT+SMALL).GT.P) GOTO 200
         BOT = CUT
      ELSE
         TOP = CUT
 
      ENDIF
      GOTO 100
 
  200 CONTINUE
      HI = CUT
 
 
 
C     Now find the bottom value satisfying....
 
      TOP = HI
      BOT = HI - 2.0
      START = .TRUE.
      IF (BOT.LT.0.) BOT = 0.
 
  300 CONTINUE
      CUT = BOT + (TOP-BOT)/2
 
      IF (START) THEN
         START = .FALSE.
      ELSEIF (CUT.EQ.OLD) THEN
         GOTO 400
      ENDIF
 
      OLD = CUT
 
      IF (CUT.EQ.HI) THEN
         X = HI
         IER = 0
         RETURN
      ENDIF
 
      ABMERFC = ABMERF(CUT)
 
 
 
      IF (ABMERFC.GE.P) THEN
         IF (ABMERF(CUT-SMALL).LT.P) GOTO 400
         TOP = CUT
      ELSE
         BOT = CUT
 
      ENDIF
      GOTO 300
 
 
  400 CONTINUE
      X = CUT + (HI-CUT)/2
 
      IER = 0

      RETURN
      END

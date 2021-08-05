CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GSER(GAMSER,A,X,GLN)
C                ********************
C     Returns the incomplete gamma function P(a,x) as series GAMSER.
C     Also returns ln gamma(a) as GLN
 
      INTEGER ITMAX, N
      REAL EPS,SUM,X,GAMSER,GLN
      REAL A,AP,GAMMLN,DEL
 
      PARAMETER (ITMAX = 10000, EPS = 3.E-20)
 
COM is this a function call or what?
      GLN = GAMMLN(A)
 
      IF (X.LE.0.0) THEN
         IF (X.LT.0.0) THEN
            PRINT *, 'Error in GSER'
            CALL DIE
         ENDIF
         GAMSER = 0.0
         RETURN
      ENDIF
      AP = A
      SUM = 1.0/A
      DEL = SUM
      DO 100 N = 1, ITMAX
         AP = AP + 1.0
         DEL = DEL*X/AP
         SUM = SUM + DEL
         IF (ABS(DEL).LT.ABS(SUM)*EPS) GOTO 200
  100 CONTINUE
      PRINT *, 'Error in GSER: A too large, ITMAX too small'
      CALL DIE
  200 CONTINUE
      GAMSER = SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END

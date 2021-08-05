      FUNCTION GAMMQ(A,X)
 
COM OMUPD BNJ 2/9/91
 
      REAL X,A,GLN,GAMMQ,GAMSER,GAMMCF
 
COM
 
C     Returns the incomplete gamma function Q(a,p) == 1 - P(a,x)
 
      IF (X.LT.0..OR.A.LE.0.) THEN
         PRINT *, 'Error in GAMMQ'
         CALL DIE
      ENDIF
      IF (X.LT.A+1.) THEN
         CALL GSER(GAMSER,A,X,GLN)
         GAMMQ = 1. - GAMSER
      ELSE
         CALL GCF(GAMMCF,A,X,GLN)
         GAMMQ = GAMMCF
      ENDIF
      RETURN
      END

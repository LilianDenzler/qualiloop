CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      FUNCTION ABMERFC(X)
C              *******
C     Returns the complementary error function erfc(x)
C
 
COM OMUPD BNJ
      REAL X,ABMERFC,GAMMP,GAMMQ
COM
      IF (X.LT.0.) THEN
         ABMERFC = 1. + GAMMP(.5,X**2)
      ELSE
         ABMERFC = GAMMQ(.5,X**2)
      ENDIF
      RETURN
      END

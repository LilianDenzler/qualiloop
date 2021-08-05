      REAL FUNCTION INTERPA(X1,Y1,X2,Y2,X)
C
C     Performs linear interpolation between the points X1,Y1 and X2,Y2
C     for X. The value of the INTERPOLATE is Y. The interpolation
C     formula attempts to postpone until the last moment the subtraction
C     of roughly equal magnitude quantities.
C
C      IMPLICIT INTEGER(A-Z)
      IMPLICIT NONE
 
      REAL X1, Y1, X2, Y2, X
      REAL*8 M
C
      IF (X.EQ.X1) THEN
         INTERPA = Y1
      ELSEIF (X.EQ.Y2) THEN
         INTERPA = Y2
      ELSE
         M = (DBLE(Y2)-Y1)/(DBLE(X2)-X1)
         INTERPA = M*(DBLE(X)-X1) + Y1
      ENDIF
 
      RETURN
      END

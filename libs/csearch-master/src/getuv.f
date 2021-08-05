      SUBROUTINE GETUV(X,Y,Z,IND,U,V)
C      IMPLICIT INTEGER(A-Z)
      IMPLICIT NONE
 
      REAL X(*), Y(*), Z(*), U(*), V(*)
      REAL DX(3), R, DOT, VLEN, DOTUR
      INTEGER IND(3),I
 
      DX(1) = X(IND(2)) - X(IND(1))
      DX(2) = Y(IND(2)) - Y(IND(1))
      DX(3) = Z(IND(2)) - Z(IND(1))
      R = VLEN(DX,3)
      CALL ASCALE(DX,1.0/R,U,3)
      DX(1) = X(IND(3)) - X(IND(2))
      DX(2) = Y(IND(3)) - Y(IND(2))
      DX(3) = Z(IND(3)) - Z(IND(2))
      DOTUR = DOT(U,DX,3)
      DO 100 I = 1, 3
         V(I) = DX(I) - DOTUR*U(I)
  100 CONTINUE
      R = VLEN(V,3)
      CALL ASCALE(V,1.0/R,V,3)
      RETURN
      END

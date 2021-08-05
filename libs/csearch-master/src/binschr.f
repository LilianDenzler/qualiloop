CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BINSCHR(A,N,KEY,ISTART,ISTOP)
C
C     Performs a binary search on the array A for the the KEY. ISTART and
C     ISTOP will return the array elements which bound KEY. ISTART will be
C     zero if KEY precedes the first element; ISTOP will be greater than
C     N if KEY comes after the last element.
 
      IMPLICIT NONE

C -- Temp
      include "params.inc"
      include "coords.inc"
      include "getpar.inc"
      include "sidetop.inc"
C -- End Temp


COM OMUPD BNJ 29/08/91
      INTEGER N,ISTART,ISTOP,IND
COM
 
      REAL A(*), KEY
C
      IF (N.EQ.0) THEN
         ISTART = 0
         ISTOP = 1
      ELSEIF (KEY.LT.A(1)) THEN
         ISTART = 0
         ISTOP = 1
      ELSEIF (KEY.GT.A(N)) THEN
         ISTART = N
         ISTOP = N + 1
      ELSE
         ISTART = 1
         ISTOP = N
   50    CONTINUE
         IF (ISTOP-ISTART.GT.1) THEN
            IND = (ISTOP+ISTART)/2
            IF (KEY.LE.A(IND)) THEN
               ISTOP = IND
            ELSE
               ISTART = IND
            ENDIF
            GOTO 50
         ENDIF
      ENDIF
      RETURN
      END
 

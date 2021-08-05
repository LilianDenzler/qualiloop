CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INFRLS(HEAD,TAIL,NEXT,N)
C
C     Initializes an array into a free list, i.e. everything in the
C     array is linked sequentially into one long list. The length of the
C     array is given by N.
C
C      IMPLICIT INTEGER(A-Z)
      IMPLICIT NONE
 
COM OMUPD BNJ 29/08/91
      INTEGER HEAD,N,TAIL,I
COM
 
      INTEGER NEXT(*)
C
      IF (N.EQ.0) THEN
         HEAD = 0
         TAIL = 0
      ELSE
         HEAD = 1
         TAIL = N
         DO 50 I = 2, N
            NEXT(I-1) = I
   50    CONTINUE
         NEXT(N) = 0
      ENDIF
      RETURN
      END

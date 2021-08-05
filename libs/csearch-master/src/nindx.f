CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NINDX
C     Recoded ACRM 12.06.91
C     Finds NUMBER in sorted NARRAY (length NLEN) by binary search
C     Returns its index in the array, 0 if not found.
 
      INTEGER FUNCTION NINDX(NUMBER,NARRAY,NLEN)
      INTEGER NUMBER, NARRAY(*), NLEN
      INTEGER ISTART, ISTOP
 
      ISTART = 1
      ISTOP = NLEN
 
  100 CONTINUE
      NINDX = ISTART + (ISTOP-ISTART)/2
      IF (NARRAY(NINDX).EQ.NUMBER) THEN
         RETURN
      ELSEIF (NARRAY(NINDX).GT.NUMBER) THEN
         ISTOP = NINDX - 1
      ELSEIF (NARRAY(NINDX).LT.NUMBER) THEN
         ISTART = NINDX + 1
      ENDIF
      IF (ISTOP.GE.ISTART) GOTO 100
 
      NINDX = 0
      RETURN
      END

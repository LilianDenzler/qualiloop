      SUBROUTINE STUPHB1(DONP,ACCP)
C
C     Sets up the hydrogen bond data for CGEN. DONP and ACCP are
C     arrays which mark each atom as being a donor or acceptor,
C     respectively. Non-zero elements serve as markers. Other variables
C     pertaining to hydorgen bonds are also set.
C
      IMPLICIT NONE
 
      include "params.inc"
      include "values.inc"
      include "pstruct.inc"
      include "hbonds.inc"
      include "cg.inc"
 
      INTEGER DONP(*), ACCP(*)
      INTEGER I
      CHARACTER*100 BUFFER
 
C
      CALL FILL4(DONP,NATOMS,0)
      CALL FILL4(ACCP,NATOMS,0)
 
C ACRM Do we actually need this??
      DO 100 I = 1, NDONAT-1
         IF (HBDHYD(I).EQ.0) THEN
            WRITE (BUFFER, 1000) I, NDONAT
 1000       FORMAT('Error in CGEN: HBond atom ',I5,' of ',I5,
     +             ' was undefined')
            CALL CPRINT(BUFFER)
            CALL DIE
         ENDIF
         DONP(HBDHYD(I)) = I
  100 CONTINUE
 
      DO 200 I = 1, NACCAT
         ACCP(HBACPT(I)) = I
  200 CONTINUE
      CUTHB2 = HBCUT*HBCUT
      COS_CUTHBA = COS(HBACUT*DTORAD)
 
      RETURN
      END

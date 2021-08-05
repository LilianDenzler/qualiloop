CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RNHBNDSRCH(STARTX,STARTY,STARTZ,LASTX,LASTY,
     +                            LASTZ,SPACE_GRID,CLSHD,CLSATM,
     +                            EXCLUDED,SEARCH_MODE,DONP,NEXTCLS,
     +                            PARM_NO,EHB,ACCP,IND)
 
 
      IMPLICIT NONE
      include "params.inc"
      include "grid.inc"
      include "cg.inc"
      include "coords.inc"
      include "engpar.inc"
      include "pstruct.inc"
      include "dbg.inc"
      INTEGER STARTX, STARTY, STARTZ, LASTX, LASTY, LASTZ
      INTEGER*2 SPACE_GRID(NGRIDX,NGRIDY,NGRIDZ), PARM_NO(100,100)
      INTEGER JX, JY, JZ, ISPACE, PSN, IATM, IHBOND, JHBOND, KHBOND
      INTEGER DONP(*), NEXTCLS(*), CLSHD(*), CLSATM(*), IND
      LOGICAL EXCLUDED(*)
      INTEGER SEARCH_MODE, ACCP(*)
      DOUBLE PRECISION EHB
 
      DO 100 JX = STARTX, LASTX
         DO 50 JY = STARTY, LASTY
            DO 20 JZ = STARTZ, LASTZ
               ISPACE = SPACE_GRID(JX,JY,JZ)
               IF (ISPACE.NE.0) THEN
                  PSN = CLSHD(ISPACE)
    5             CONTINUE
                  IF (PSN.NE.0) THEN
                     IATM = CLSATM(PSN)
                     IF (.NOT.EXCLUDED(IATM)) THEN
                        IF (SEARCH_MODE.EQ.DONOR_HBOND_ENERGY.AND.
     +                      ACCP(IATM).NE.0) THEN
                           IHBOND = HBDONR(DONP(IND))
                           JHBOND = IATM
                           KHBOND = IND
C                             EVALUATE-ENERGY-FOR-ONE-HYDROGEN-BOND
                           CALL HBENERGY(IHBOND,JHBOND,KHBOND,
     +                           PARM_NO,EHB)
                        ELSEIF (SEARCH_MODE.EQ.ACCEPTOR_HBOND_ENERGY
     +                          .AND.DONP(IATM).NE.0) THEN
                           IHBOND = HBDONR(DONP(IATM))
                           JHBOND = IND
                           KHBOND = IATM
C                             EVALUATE-ENERGY-FOR-ONE-HYDROGEN-BOND
                           CALL HBENERGY(IHBOND,JHBOND,KHBOND,
     +                           PARM_NO,EHB)
                        ENDIF
                     ENDIF
                     PSN = NEXTCLS(PSN)
                     GOTO 5
                  ENDIF
               ENDIF
   20       CONTINUE
   50    CONTINUE
  100 CONTINUE
      RETURN
      END

 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RNNBSRCH(SEARCH_MODE,STARTX,LASTX,STARTY,LASTY,
     +                            STARTZ,LASTZ,SPACE_GRID,CLSHD,CLSATM,
     +                            EXCLUDED,PARM_NO,MAXEVDW,IMPACT,IND,
     +                            RESBYA,QSIDE,SIDEHITS,EEL,IGNORE_EVDW,
     +                            EVDW,NEXTCLS,XIND,YIND,ZIND)
 
      IMPLICIT NONE
      include "params.inc"
      include "grid.inc"
      include "cg.inc"
      include "coords.inc"
      include "engpar.inc"
      include "pstruct.inc"
      include "dbg.inc"
      include "values.inc"
      REAL CUT2_SEARCH, CGIND
      INTEGER ITCIND, OFFIND, SEARCH_MODE, JX, JY, JZ, ISPACE
      INTEGER STARTX, STARTY, STARTZ, LASTX, LASTY, LASTZ
      INTEGER*2 SPACE_GRID(NGRIDX,NGRIDY,NGRIDZ)
      INTEGER CLSHD(*), CLSATM(*)
      LOGICAL EXCLUDED(*)
      INTEGER*2 PARM_NO(100,100)
      REAL MAXEVDW, EEL1, ENB, DELX, DELY, DELZ, XIND, YIND, ZIND
      LOGICAL IMPACT
      INTEGER IND, ITCIATM, IATM, PSN
      INTEGER RESBYA(*)
      LOGICAL QSIDE(*), IGNORE_EVDW
      INTEGER SIDEHITS(*), NEXTCLS(*)
      DOUBLE PRECISION EEL, EVDW
      REAL R2, R4, R6, R12, F1, F2, SS
      INTEGER IC, INDSEG, GETSEG, IATMSEG
      CHARACTER*100 BUFFER
 
      IF (SEARCH_MODE.EQ.CONTACT) THEN
         CUT2_SEARCH = GRID2
      ELSE
         CUT2_SEARCH = CUTNB2
      ENDIF
      CGIND = ATCHRG(IND)
      ITCIND = ATFLAG(ATCODE(IND))
      OFFIND = IOFF(ITCIND)
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
                        DELX = XCART(IATM) - XIND
                        DELY = YCART(IATM) - YIND
                        DELZ = ZCART(IATM) - ZIND
                        SS = DELX*DELX + DELY*DELY + DELZ*DELZ
                        IF (SS.LT.CUT2_SEARCH) THEN
                           IF (SEARCH_MODE.EQ.CONTACT) THEN
                              R2 = 1.0/SS
                              R4 = R2*R2
                              R6 = R4*R2
                              R12 = R6*R6
                              ITCIATM = ATFLAG(ATCODE(IATM))
                              IC = PARM_NO(ITCIND,ITCIATM)
                              F1 = VDWR12(IC)*R12
                              F2 = VDWR6(IC)*R6
                              ENB = F1 - F2
                              IF (ENB.GT.MAXEVDW) THEN
                                 IMPACT = .TRUE.
                                 IF (DBG_CGEN.GT.1) THEN
                                    INDSEG = GETSEG(RESBYA(IND),SEGNDX,
     +                                    NSEGS)
                                    IATMSEG = GETSEG(RESBYA(IATM),
     +                                    SEGNDX,NSEGS)
                                    WRITE (BUFFER,9001) SQRT(SS), ENB
                                    CALL CPRINT(BUFFER)
                                    WRITE (BUFFER,9002) SEGID(INDSEG),
     +                                    RESNME(RESBYA(IND)),
     +                                    RESID(RESBYA(IND)), 
     +                                    ATMNME(IND),
     +                                    SEGID(IATMSEG),
     +                                    RESNME(RESBYA(IATM)),
     +                                    RESID(RESBYA(IATM)),
     +                                    ATMNME(IATM)
                                    CALL CPRINT(BUFFER)
                                 ENDIF
 
                                 IF (QSIDE(IATM).AND.QSIDE(IND)) THEN
                                    IF (DBG_CGEN.GT.1) THEN
                                       WRITE (BUFFER,9003) ATMNME(IND),
     +                                       IND, ATMNME(IATM), IATM
                                       CALL CPRINT(BUFFER)
                                    ENDIF
                                    SIDEHITS(RESBYA(IATM)) = IND
                                 ENDIF
                                 GOTO 200
                              ENDIF
                           ELSE
                              R2 = 1.0/SS
                              IF (CONS_DIE) THEN
                                 EEL1 = 332.0716*CGIND*ATCHRG(IATM)
     +                                  *SQRT(R2)
     +                                  /EPSILON
                              ELSE
                                 EEL1 = 332.0716*CGIND*ATCHRG(IATM)
     +                                  *R2
                              ENDIF
                              EEL = EEL + EEL1
                              IF (.NOT.IGNORE_EVDW) THEN
                                 R4 = R2*R2
                                 R6 = R4*R2
                                 R12 = R6*R6
                                 ITCIATM = ATFLAG(ATCODE(IATM))
                                 IC = PARM_NO(ITCIND,ITCIATM)
                                 F1 = VDWR12(IC)*R12
                                 F2 = VDWR6(IC)*R6
                                 ENB = F1 - F2
                                 EVDW = EVDW + ENB
                              ENDIF
                              IF (DBG_CGEN.GT.3) THEN
                                 WRITE (BUFFER,9004) ATMNME(IND), IND,
     +                                  ATMNME(IATM), IATM, SQRT(SS),
     +                                  ENB, EEL1
                                 CALL CPRINT(BUFFER)
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                     PSN = NEXTCLS(PSN)
                     GOTO 5
                  ENDIF
               ENDIF
   20       CONTINUE
   50    CONTINUE
  100 CONTINUE
  200 CONTINUE
      RETURN
 9001 FORMAT (' Close contact failed. R =',F5.2,' E = ',1PG12.5)
 9002 FORMAT (' Atoms:',4(1X,A4),' and',4(1X,A4))
 9003 FORMAT (' Sidechain collision, atoms:',2(1X,A4,I5))
 9004 FORMAT (' Energy for atoms:',2(1X,A4,I5),' R =',F5.2,' EVDW = ',
     +        1PG12.5,' ELEC = ',1PG12.5)
      END

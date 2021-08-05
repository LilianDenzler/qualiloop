      SUBROUTINE ECCHKCONT(IATM,ATOMNO,ORI,ECX,ECY,ECZ,AVOIDX,
     +                     AVOIDR,AVOIDPHI,RCONTACT,STARTPHI,
     +                     LASTPHI,IMPACT,SIDEHITS,RESBYA,CNTNBX,
     +                     NBXA,QSIDE)
C
C     Determines the range of torsion angles for the free torsion in the
C     clump of which atom number, ATOMNO, is a member, where there is no
C     contact closer than those found in RCONTACT to the atom numbered
C     IATM. ORI, ECX, ECY, and ECZ define the local coordinate system
C     for the clump. AVOIDX gives the distance along the x axis for
C     ATOMNO, AVOIDR gives the rotation radius for that atom, AVOIDPHI
C     gives the offset from the free torsion, RCONTACT is an array of
C     van der Waals radii for all atom types with ATOMNO. STARTPHI and
C     LASTPHI give the allowed range where STARTPHI is always less than
C     or equal to LASTPHI and STARTPHI is within [-pi,pi]. If no
C     placement is possible, then STARTPHI is returned with -10.0. If
C     all placements are possible, then IMPACT is turned off. The reason
C     for this redundancy over STARTPHI and LASTPHI is the imprecision
C     of floating point numbers. SIDEHITS records any sidechain hits for
C     ATOMNO in the event that both IATM and ATOMNO are sidechain atoms
C     and IATM could make close contact. Finally, RESBYA gives residue
C     numbers by atom to permit quick access into SIDEHITS.
C
 
C      IMPLICIT INTEGER(A-Z)
 
COM OMUPD BNJ 2/9/91
 
      IMPLICIT NONE
 
COM
 
      include "params.inc"
      include "cg.inc"
      include "values.inc"
      include "coords.inc"
      include "dbg.inc"
      include "grid.inc"
      include "engpar.inc"
      include "pstruct.inc"
C
      INTEGER I1,I
 
      INTEGER IATM, ATOMNO
      REAL ORI(3), ECX(3), ECY(3), ECZ(3), AVOIDX, AVOIDR, AVOIDPHI
      REAL RCONTACT(*), STARTPHI, LASTPHI
      LOGICAL IMPACT
      INTEGER SIDEHITS(*), RESBYA(*), CNTNBX(*)
      INTEGER*2 NBXA(MAXNBX,*)
      LOGICAL QSIDE(*)
C
      REAL XT, YT, ZT, XS, YS, ZS, DELTA1, DELTAP, AA, BB, C, D, CC, DD
COM   ACRM Change BYTE to CHARACTER
COM      BYTE ST1(21),ST2(21)
      CHARACTER*21 ST1, ST2
 
      CHARACTER*100 BUFFER
C
      STARTPHI = -PI
      LASTPHI = PI
      IMPACT = .FALSE.
      DO 100 I = 1, CNTNBX(ATOMNO)
         IF (NBXA(I,ATOMNO).EQ.IATM) RETURN
  100 CONTINUE
C
C     Get coordinates of IATM in coordinate system for ATOMNO.
C
      XT = XCART(IATM) - ORI(1)
      YT = YCART(IATM) - ORI(2)
      ZT = ZCART(IATM) - ORI(3)
      XS = XT*ECX(1) + YT*ECX(2) + ZT*ECX(3)
      YS = XT*ECY(1) + YT*ECY(2) + ZT*ECY(3)
      ZS = XT*ECZ(1) + YT*ECZ(2) + ZT*ECZ(3)
      I1 = ATFLAG(ATCODE(IATM))
      AA = RCONTACT(I1)**2 - (AVOIDX-XS)**2
      BB = AA - AVOIDR**2 - YS**2 - ZS**2
      C = -BB/2.0/AVOIDR
      CC = C*C
      DD = YS**2 + ZS**2
      IF (DD.LT.CC.AND.C.GT.0.0) THEN
      ELSEIF (DD.LE.CC.AND.C.LE.0.0) THEN
         STARTPHI = -10.0
         IMPACT = .TRUE.
C        CHECK-SIDEHITS
         IF (QSIDE(IATM).AND.QSIDE(ATOMNO)) THEN
            IF (DBG_CGEN.GT.1) THEN
               CALL CPRINT(' Next contact entered in sidehits array.')
            ENDIF
            SIDEHITS(RESBYA(IATM)) = ATOMNO
         ENDIF
      ELSE
         IMPACT = .TRUE.
         D = SQRT(DD)
         DELTA1 = ACOS(C/D)
         DELTAP = ATAN2(ZS,YS) - AVOIDPHI
         STARTPHI = DELTAP + DELTA1
         LASTPHI = DELTAP - DELTA1 + 2*PI
  150    CONTINUE
         IF (STARTPHI.GT.PI) THEN
            STARTPHI = STARTPHI - 2*PI
            LASTPHI = LASTPHI - 2*PI
            GOTO 150
         ENDIF
  200    CONTINUE
         IF (STARTPHI.LT.-PI) THEN
            STARTPHI = STARTPHI + 2*PI
            LASTPHI = LASTPHI + 2*PI
            GOTO 200
         ENDIF
C        CHECK-SIDEHITS
         IF (QSIDE(IATM).AND.QSIDE(ATOMNO)) THEN
            IF (DBG_CGEN.GT.1) THEN
               CALL CPRINT(' Next contact entered in sidehits array.')
            ENDIF
            SIDEHITS(RESBYA(IATM)) = ATOMNO
         ENDIF
      ENDIF
      IF (DBG_CGEN.GT.1) THEN
         IF (DBG_CGEN.GT.2.OR.STARTPHI.NE.-PI.OR.LASTPHI.NE.PI) THEN
 
COM ACRM Changed to use FORPRTATM rather then BYTES...
            CALL FORPRTATM(ST1,IATM)
            CALL FORPRTATM(ST2,ATOMNO)
 
            WRITE (BUFFER,9001) ST1, ST2
            CALL CPRINT(BUFFER)
 
            WRITE (BUFFER,9002) STARTPHI/DTORAD, LASTPHI/DTORAD
            CALL CPRINT(BUFFER)
         ENDIF
         IF (DBG_CGEN.GT.2) THEN
            WRITE(BUFFER,9003) ORI,XCART(IATM),YCART(IATM),ZCART(IATM)
            CALL CPRINT(BUFFER)
         ENDIF
      ENDIF
      RETURN
 9001 FORMAT (' EC_CHECK_CONTACT: ',A20,' and ',A20)
 9002 FORMAT ('           Angles: ',2F10.3,' degrees.')
 9003 FORMAT (' Origin is at ',3F8.3,'. Clump atom: ',3F8.3,'.')
      END

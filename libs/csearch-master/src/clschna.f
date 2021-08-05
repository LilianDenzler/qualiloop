COM OMUPD BNJ 19/11/91
COM
COM CISTRANS ==> ENTZUS
COM PEPOMEGA ==> POMEGA
COM NEWANGLE ==> NANGLE
COM ITERMAX ==> ITMAX
COM ITERPHI ==> ITRPHI
COM DIFFPHI ==> DIFPHI
COM GLOWPHI ==> GLOPHI
COM BASEMAXIT ==> BSMXIT
COM VARDEBUG ==> VARDBG
COM MAXITRF ==> MXITRF
COM OLD_NP ==> NPOLD
COM SMALLEST_NP ==> NPLITL
COM REFER_NP ==> REFNP
COM NEW_NP ==> NPNEW
COM SMALLEST_SGNS ==> SGNSML
COM REFER_SGNS ==> SGNREF
COM NPOLDHI1 ==> NPOLH1
COM SAVED_I* ==> SAVI*
COM MODIFY_TRIES ==> MODTRI
COM NSDSTEP ==> NSDSTP
COM MAXITMIN ==> MXITMN
COM MAX_MODIFY_TRIES ==> MXMODT
COM BIG_G ==> BIGGEE
COM SMALLEST_G ==> SMLGEE
COM SMALLEST_OMEGA1 ==> SMLOM1
COM OLD_PHI1 ==> OLDP1
COM PHIRANGE ==> PHIRNG
COM PHI1_FRAC ==> P1FRAC
COM REFER_OMEGA1 ==> REFOM1
COM DISPLAY_SG ==> DISPSG
COM GPHIPREV ==> GPPREV
COM SAVED_ANGLE ==> SAVANG
COM DOTMOVEGRAD ==> DTMVGD
COM NORMSHIFT ==> NMSHFT
COM LENSHIFT ==> LNSHFT
COM PREVANGLE ==> PRVANG
COM LENMOVE ==> LENMV
COM POS_SEEN ==> POSSEN
COM NEG_SEEN ==> NEGSEN
COM ZERO_G_FOUND ==> FND0G
COM GBOUND_OK ==> GBNDOK
COM MODIFY_DONE ==> MODDUN
COM OLDANGLE ==> OLDANG
COM OLD_OMEGA1 ==> OLDOM1
COM STOP_SD ==> STOPSD
COM DROP_STEP ==> DRPSTP
COM GO_ON ==> CARYON
COM ADJUSTED_THETA ==> ADTHET
COM REF1PHI ==> RF1PHI
COM REF2PHI ==> RF2PHI
COM NREF1PHI ==> NRF1PY
COM NREF2PHI ==> NRF2PY
COM DGDANGLE ==> DGDANG
COM IRESTMAX ==> IRSTMX
COM
COM
      SUBROUTINE CLSCHNA(X,Y,Z,ATMIND,BOND,ANGLE,NEWX,NEWY,NEWZ,OMEGA,
     +                  ITER,ENTZUS,POMEGA,MAXDT,NANGLE,MAXG)
COM
COM      IMPLICIT CHARACTER*1000(A-H,J-Z)
      IMPLICIT NONE
COM
C
C
      INTEGER ITMAX
      INTEGER ITER
      INTEGER ATMIND(3,7)
      INTEGER ITERCT, ITRPHI
      INTEGER SGNS(2), NROOT, NP
C
      INTEGER MAXRTS
      PARAMETER (MAXRTS=50)
C
COM OMLUPD BNJ
C
      INTEGER I99929, I99913  ,I99909  ,I99862  
      INTEGER IRSTMX, IREST   ,I99856  ,I99851  ,I99848  ,I99818  
      INTEGER I99821, I99815  ,I99811  ,I99803  ,I99805  ,I99780  
      INTEGER I99765, I99703  ,I99705  ,I99699  ,I99678  ,IER     
      INTEGER I99668, I99666  ,I99662  ,I99660  ,I99632  ,I99630  
C
      INTEGER I99998, I99996, I99994, I     , I99989, I99987
      INTEGER I99985, I99981, I99979, I99977, IS2   , IS1
      INTEGER I99969, I99975, ILIM  , I99956, I99628, I99626
      INTEGER I99624, I99612, I99610, I99608, I99602, NINT,IND
C
      INTEGER KOUNT, IBEG, IEN, KN, N, J, MAXIT, BSMXIT
      INTEGER DEBUG, VARDBG
      INTEGER MXITRF, MAXNW
      INTEGER NPOLD, NPLITL, REFNP, NPNEW
      INTEGER SGNSML(2)
      INTEGER SGNREF(2)
      INTEGER NPOLH1, IT1, IT2, IP
      INTEGER SAVI1, SAVI2, SAVI3, SAVI4, SAVI5, SAVI6
      INTEGER SAVI7, SAVI8, SAVI9
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I9
      INTEGER MODTRI, NSDSTP
      INTEGER MXITMN, MXMODT
      INTEGER NPHI1
C
C
      REAL X(*), Y(*), Z(*)
      REAL BOND(3,3), ANGLE(3,3), NEWX(6), NEWY(6), NEWZ(6), OMEGA(6)
      REAL S(3), U(3), V(3), DX(3), U0(3), V0(3), W0(3)
      REAL U6(3), V6(3), Q(3,0:2), R(3), RHO(0:2), SIGMA(0:2)
      REAL THETA(0:5), CTHETA(0:5), STHETA(0:5)
      REAL P(3,0:5)
      REAL POMEGA(3), MAXDT, NANGLE(3,3), MAXG
C
C
      REAL BIGGEE
      PARAMETER (BIGGEE=10.0)
C
C
      REAL PHIRTS(MAXRTS), RTS(MAXRTS)
      REAL LOWPHI, HIPHI, DIFPHI, GLOPHI, GHIPHI
      REAL RT, H, DELFPR, FRTDEF, LAMBDA, DELF, DFPRLM, NUM
      REAL DEN, G, SQR, FRT, FRTPRV, FRTDEFPREV
      REAL FRTDEFPOS, FRTDEFNEG, OMEGA1POS, OMEGA1NEG
      REAL OMPRV1, OMPRV2
      REAL ABMERF, FRF, GRF, WRF, SFRF, SWRF, ARF, BRF, RATIO
      REAL FX, DFX, FXDFX, RF1PHI1, RF2PHI1, NRF1PY1, NRF2PY1
      REAL DW, DXX, TREF
      REAL DGDANG(3,3), OLDANG(3,3)
      REAL DT, SG1, SG2, SMLGEE, ANG1, ANG2, ABSG, SMLOM1
      REAL SG3, OLDP1(2,2)
      REAL FRAC1, FRAC2, PHIRNG, DIFF1, DIFF2, P1FRAC
      REAL OLDOM1, REFOM1, GPPREV
      REAL SAVANG, DISPSG
      REAL OLDG, DTMVGD, LENSD, STEPSD, NMSHFT(3,3), LNSHFT
      REAL PRVANG(3,3), LENMV
      REAL TA(3), TB(3), SOMEGA(6), COMEGA(6)
      REAL PHI1(2,2), GPHI1(2,2,-1:1,-1:1)
      REAL RR2, RR, DOT, T1, T2, T3, T4, DENOM, F1R, F2R, F3R, F4R, UB,
     +     LB
      REAL DISC1, DISC2
      REAL QR12, QR22, QR1, QR2, RQSUM
      REAL TEMPA(3), T(3), A, B, C, DELTA
      REAL APHI, BPHI
      REAL GPHI, PRODG
      REAL XX, YY, ZZ, W, Q3(3)
      REAL G1, G2, G3, G4
      REAL DIFF12, DIFF13, DIFF14
      REAL DIFF23, DIFF24, DIFF34
      REAL EPS1, EPS2
      PARAMETER (EPS1=1.0E-6)
      PARAMETER (EPS2=1.0E-6)
      REAL EPS3
      REAL RTEMP1,RTEMP2
C
C
      CHARACTER BUFFER*200
C
C
      LOGICAL ENTZUS(3)
      LOGICAL OVERLP, POSSEN, NEGSEN, MULLER
      LOGICAL DONERF
      LOGICAL DONENW, TOOBIG
      LOGICAL FND0G, GBNDOK
      LOGICAL MODDUN, ABORT, RETRY_CLSCHN, OK
      LOGICAL ADTHET
      LOGICAL STOPSD, DRPSTP
      LOGICAL SWAP34
      LOGICAL FOUND, DONE, CARYON
C
C
C
      INCLUDE "values.inc"
      INCLUDE "dbg.inc"
C
C
C
      DATA RTS/MAXRTS*0.0/
      DATA MXITRF/50/
      DATA MAXNW/5/
      DATA MXMODT/4/
      DATA MXITMN/4/
      DATA BSMXIT/10/
      DATA EPS3/2.0E-3/
      DATA ITMAX/0/
C
C
C
      DEBUG = MOD(DBG_CLSCHN,10)
      VARDBG = MOD(DBG_CLSCHN/10,10)
C
      IF (ITER.EQ.0) THEN
         ASSIGN 100 TO I99998
         ASSIGN 300 TO I99994
         GOTO 2700
      ENDIF
  100 CONTINUE
      ASSIGN 200 TO I99996
      GOTO 4200
  200 CONTINUE
      RETURN
  300 CONTINUE
      DO 400 I = 1, 3
         IF (.NOT.(ENTZUS(I))) THEN
            POMEGA(I) = PI
         ELSE
            POMEGA(I) = 0.0
         ENDIF
  400 CONTINUE
      ASSIGN 500 TO I99989
      GOTO 800
  500 CONTINUE
      ASSIGN 600 TO I99987
      GOTO 1000
  600 CONTINUE
      ADTHET = .FALSE.
      SMLGEE = BIGGEE
      ASSIGN 700 TO I99985
      GOTO 1600
  700 CONTINUE
      SGNS(1) = 1
      SGNS(2) = 1
      NROOT = 0
      NP = 1
      ITERCT = 0
      ITRPHI = 0
      GOTO I99998
  800 CONTINUE
      DO 900 I1 = 1, 3
         DO 850 I2 = 1, 3
            NANGLE(I1,I2) = ANGLE(I1,I2)
  850    CONTINUE
  900 CONTINUE
      GOTO I99989
 1000 CONTINUE
      ASSIGN 1100 TO I99981
      GOTO 2100
 1100 CONTINUE
      ASSIGN 1200 TO I99979
      GOTO 2300
 1200 CONTINUE
      ASSIGN 1300 TO I99977
      GOTO 2800
 1300 CONTINUE
      GOTO I99987
 1400 CONTINUE
      DISPSG = SMLGEE
      SMLGEE = 0.0
      DO 1500 NP = 1, NPHI1
         DO 1450 IS2 = 1, 2
            DO 1440 IS1 = 1, 2
               SGNS(1) = 3 - 2*IS1
               SGNS(2) = 3 - 2*IS2
               OMEGA(1) = PHI1(1,NP)
 1410          CONTINUE
               IF (OMEGA(1).GT.PHI1(2,NP)) THEN
                  OMEGA(1) = PHI1(2,NP)
                  ASSIGN 1430 TO I99969
               ELSE
                  ASSIGN 1420 TO I99969
               ENDIF
               GOTO 14100
 1420          CONTINUE
               WRITE (BUFFER,'(2PG14.7)') OMEGA(1)/DTORAD, GPHI
               CALL CPRINT(BUFFER)
               OMEGA(1) = OMEGA(1) + .02
               GOTO 1410
 1430          CONTINUE
               WRITE (BUFFER,'(2PG14.7)') OMEGA(1)/DTORAD, GPHI
               CALL CPRINT(BUFFER)
 1440       CONTINUE
 1450    CONTINUE
 1500 CONTINUE
      SMLGEE = DISPSG
      GOTO I99975
 1600 CONTINUE
      DO 1800 NP = 1, NPHI1
         DO 1650 IS1 = 1, 2
            SGNS(1) = 3 - 2*IS1
            DO 1620 IS2 = 1, 2
               SGNS(2) = 3 - 2*IS2
               DO 1610 ILIM = 1, 2
                  OMEGA(1) = PHI1(ILIM,NP)
                  ASSIGN 1605 TO I99969
                  GOTO 14100
 1605             CONTINUE
                  GPHI1(ILIM,NP,SGNS(1),SGNS(2)) = GPHI
 1610          CONTINUE
 1620       CONTINUE
 1650    CONTINUE
         IF (PHI1(1,NP).NE.-PI.OR.PHI1(2,NP).NE.PI) THEN
            DO 1700 ILIM = 1, 2
               G1 = GPHI1(ILIM,NP,1,1)
               G2 = GPHI1(ILIM,NP,-1,1)
               G3 = GPHI1(ILIM,NP,1,-1)
               G4 = GPHI1(ILIM,NP,-1,-1)
               DIFF12 = G1 - G2
               DIFF13 = G1 - G3
               DIFF14 = G1 - G4
               DIFF23 = G2 - G3
               DIFF24 = G2 - G4
               DIFF34 = G3 - G4
               IF (ABS(DIFF12).GT.EPS3) THEN
                  IF (ABS(DIFF13).LE.EPS3) THEN
                     IF (ABS(DIFF24).LE.EPS3) GOTO 1670
                     ASSIGN 1670 TO I99956
                     GOTO 1900
                  ELSEIF (ABS(DIFF14).GT.EPS3) THEN
                     ASSIGN 1700 TO I99956
                     GOTO 1900
                  ELSE
                     IF (ABS(DIFF23).LE.EPS3) GOTO 1680
                     ASSIGN 1680 TO I99956
                     GOTO 1900
                  ENDIF
               ELSEIF (ABS(DIFF34).GT.EPS3) THEN
                  ASSIGN 1660 TO I99956
                  GOTO 1900
               ENDIF
 1660          CONTINUE
               IF (G1*G2.LT.0.0) THEN
                  IF (G1.GE.0.0) THEN
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G1 and G2 are being adjusted down.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,1) = 2*EPS2
                     GPHI1(ILIM,NP,-1,1) = 0.5*EPS2
                  ELSE
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G1 and G2 are being adjusted up.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,1) = 0.5*EPS2
                     GPHI1(ILIM,NP,-1,1) = 2*EPS2
                  ENDIF
               ENDIF
               IF (G3*G4.LT.0.0) THEN
                  IF (G3.GE.0.0) THEN
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G3 and G4 are being adjusted down.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,-1) = 2*EPS2
                     GPHI1(ILIM,NP,-1,-1) = 0.5*EPS2
                  ELSE
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G3 and G4 are being adjusted up.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,-1) = 0.5*EPS2
                     GPHI1(ILIM,NP,-1,-1) = 2*EPS2
                  ENDIF
               ENDIF
               GOTO 1700
 1670          CONTINUE
               IF (G1*G3.LT.0.0) THEN
                  IF (G1.GE.0.0) THEN
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G1 and G3 are being adjusted down.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,1) = 2*EPS2
                     GPHI1(ILIM,NP,1,-1) = 0.5*EPS2
                  ELSE
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G1 and G3 are being adjusted up.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,1) = 0.5*EPS2
                     GPHI1(ILIM,NP,1,-1) = 2*EPS2
                  ENDIF
               ENDIF
               IF (G2*G4.LT.0.0) THEN
                  IF (G2.GE.0.0) THEN
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G2 and G4 are being adjusted down.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,-1,1) = 2*EPS2
                     GPHI1(ILIM,NP,-1,-1) = 0.5*EPS2
                  ELSE
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G2 and G4 are being adjusted up.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,-1,1) = 0.5*EPS2
                     GPHI1(ILIM,NP,-1,-1) = 2*EPS2
                  ENDIF
               ENDIF
               GOTO 1700
 1680          CONTINUE
               IF (G1*G4.LT.0.0) THEN
                  IF (G1.GE.0.0) THEN
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G1 and G4 are being adjusted down.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,1) = 2*EPS2
                     GPHI1(ILIM,NP,-1,-1) = 0.5*EPS2
                  ELSE
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G1 and G4 are being adjusted up.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,1) = 0.5*EPS2
                     GPHI1(ILIM,NP,-1,-1) = 2*EPS2
                  ENDIF
               ENDIF
               IF (G3*G2.LT.0.0) THEN
                  IF (G3.GE.0.0) THEN
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G2 and G3 are being adjusted down.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,-1) = 2*EPS2
                     GPHI1(ILIM,NP,-1,1) = 0.5*EPS2
                  ELSE
                     IF (DEBUG.GT.0) THEN
                        WRITE (BUFFER,*)
     +                          'G2 and G3 are being adjusted up.'
                        CALL CPRINT(BUFFER)
                     ENDIF
                     GPHI1(ILIM,NP,1,-1) = 0.5*EPS2
                     GPHI1(ILIM,NP,-1,1) = 2*EPS2
                  ENDIF
               ENDIF
 1700       CONTINUE
         ENDIF
 1800 CONTINUE
      GOTO I99985
 1900 CONTINUE
      WRITE (BUFFER,9001) ILIM, NP, G1, G2, G3, G4
      CALL CPRINT(BUFFER)
      ASSIGN 2000 TO I99929
      GOTO 3500
 2000 CONTINUE
      GOTO I99956
 2100 CONTINUE
      DO 2200 I = 1, 3
         P(1,2*I-2) = BOND(1,I) - BOND(2,I)*COS(NANGLE(1,I))
         P(2,2*I-2) = BOND(2,I)*SIN(NANGLE(1,I))
         P(3,2*I-2) = 0.0
         P(1,2*I-1) = BOND(3,I)
         P(2,2*I-1) = 0.0
         P(3,2*I-1) = 0.0
 2200 CONTINUE
      GOTO I99981
 2300 CONTINUE
      DO 2400 I = 1, 3
         THETA(2*I-1) = PI - NANGLE(3,I)
         IF (POMEGA(I).NE.0.0) THEN
            THETA(2*I-2) = NANGLE(2,I) - NANGLE(1,I)
         ELSE
            THETA(2*I-2) = 2*PI - NANGLE(1,I) - NANGLE(2,I)
         ENDIF
 2400 CONTINUE
      DO 2500 I = 0, 5
         CTHETA(I) = COS(THETA(I))
         STHETA(I) = SIN(THETA(I))
 2500 CONTINUE
      DO 2600 I = 0, 2
         Q(1,I) = P(1,2*I+1)*CTHETA(2*I) + P(1,2*I)
         Q(2,I) = P(1,2*I+1)*STHETA(2*I) + P(2,2*I)
         Q(3,I) = 0.0
         RHO(I) = Q(1,I)
         SIGMA(I) = Q(2,I)
 2600 CONTINUE
      GOTO I99979
 2700 CONTINUE
      CALL GETUV(X,Y,Z,ATMIND(1,1),U0,V0)
      W0(1) = U0(2)*V0(3) - U0(3)*V0(2)
      W0(2) = U0(3)*V0(1) - U0(1)*V0(3)
      W0(3) = U0(1)*V0(2) - U0(2)*V0(1)
      CALL GETUV(X,Y,Z,ATMIND(1,7),U6,V6)
      DX(1) = X(ATMIND(1,7)) - X(ATMIND(1,1))
      DX(2) = Y(ATMIND(1,7)) - Y(ATMIND(1,1))
      DX(3) = Z(ATMIND(1,7)) - Z(ATMIND(1,1))
      S(1) = DOT(DX,U0,3)
      S(2) = DOT(DX,V0,3)
      S(3) = DOT(DX,W0,3)
      U(1) = DOT(U6,U0,3)
      U(2) = DOT(U6,V0,3)
      U(3) = DOT(U6,W0,3)
      V(1) = DOT(V6,U0,3)
      V(2) = DOT(V6,V0,3)
      V(3) = DOT(V6,W0,3)
      GOTO I99994
 2800 CONTINUE
      RR2 = (S(1)-Q(1,0))**2 + (S(2)-Q(2,0))**2 + (S(3)-Q(3,0))**2
      RR = SQRT(RR2)
      QR12 = Q(1,1)**2 + Q(2,1)**2 + Q(3,1)**2
      QR1 = SQRT(QR12)
      QR22 = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
      QR2 = SQRT(QR22)
      IF (DEBUG.GT.0.OR.VARDBG.GT.0) THEN
         WRITE (BUFFER,9002) RR, QR1, QR2
         CALL CPRINT(BUFFER)
      ENDIF
      DISC1 = (QR1+QR2)**2 - RR2
      DISC2 = RR2 - (QR1-QR2)**2
      IF (DISC1.GE.0.0.AND.DISC2.GE.0.0) THEN
         RQSUM = RR2 + QR12 - QR22
         DENOM = 2*QR12
         T1 = RHO(1)*RQSUM
         T2 = SIGMA(1)*SQRT(DISC1)*SQRT(DISC2)
         F1R = (T1+T2)/DENOM
         F2R = (T1-T2)/DENOM
         SWAP34 = .FALSE.
         T1 = SIGMA(2)*STHETA(3)
         IF (T1.LT.0) SWAP34 = .NOT.SWAP34
         T2 = RHO(2)*CTHETA(3)
         T3 = -RHO(1)*CTHETA(2) + RQSUM/2/SIGMA(1)*STHETA(2) - SIGMA(1)
     +        *STHETA(2)
         DENOM = RHO(1)*STHETA(2)/SIGMA(1) - CTHETA(2)
         IF (DENOM.LT.0) SWAP34 = .NOT.SWAP34
         F3R = (T3+T1-T2)/DENOM
         F4R = (T3-T1-T2)/DENOM
         IF (SWAP34) THEN
            T4 = F3R
            F3R = F4R
            F4R = T4
         ENDIF
         UB = AMIN1(F1R,F3R)
         LB = AMAX1(F2R,F4R)
         IF (DEBUG.GT.0.OR.VARDBG.GT.0) THEN
            WRITE (BUFFER,9003) LB, UB
            CALL CPRINT(BUFFER)
         ENDIF
         IF (LB.LE.UB) THEN
            DO 2820 I = 1, 3
               TEMPA(I) = S(I) - Q(I,0)
 2820       CONTINUE
            T(1) = CTHETA(0)*TEMPA(1) + STHETA(0)*TEMPA(2)
            T(2) = -STHETA(0)*TEMPA(1) + CTHETA(0)*TEMPA(2)
            T(3) = TEMPA(3)
            C = SQRT(T(2)**2+T(3)**2)
            APHI = (UB-T(1)*CTHETA(1))/(STHETA(1)*C)
            BPHI = (LB-T(1)*CTHETA(1))/(STHETA(1)*C)
            DELTA = ATAN2(T(2),T(3))
            ASSIGN 2900 TO I99913
            GOTO 3700
         ELSE
            NPHI1 = 0
            GOTO 3200
         ENDIF
      ELSE
         NPHI1 = 0
         GOTO 3300
      ENDIF
 2900 CONTINUE
      DO 3100 NP = 1, NPHI1
         IF (PHI1(1,NP).EQ.-PI.AND.PHI1(2,NP).EQ.PI) GOTO 3100
         PHI1(1,NP) = PHI1(1,NP) - DELTA
         PHI1(2,NP) = PHI1(2,NP) - DELTA
         ILIM = 1
         ASSIGN 2950 TO I99909
         GOTO 3900
 2950    CONTINUE
         IF (RF2PHI1.GE.PHI1(1,NP)) THEN
            PHI1(1,NP) = MIN(PHI1(1,NP),RF1PHI1)
         ELSE
            PHI1(1,NP) = RF2PHI1
         ENDIF
         ILIM = 2
         ASSIGN 3000 TO I99909
         GOTO 3900
 3000    CONTINUE
         IF (RF1PHI1.LE.PHI1(2,NP)) THEN
            PHI1(2,NP) = MAX(PHI1(2,NP),RF2PHI1)
         ELSE
            PHI1(2,NP) = RF1PHI1
         ENDIF
 3100 CONTINUE
      IF (DEBUG.GT.0) THEN
         WRITE (BUFFER,9004) (PHI1(1,NP),PHI1(2,NP),NP=1,NPHI1)
         CALL CPRINT(BUFFER)
      ENDIF
      IF (NPHI1.EQ.2) THEN
         IF (PHI1(2,1).GT.PHI1(1,2)) THEN
            IF (DEBUG.GT.0) THEN
               WRITE (BUFFER,'(A)') ' PHI1 overlap corrected.'
               CALL CPRINT(BUFFER)
            ENDIF
            PHI1(2,1) = PHI1(2,2)
            NPHI1 = 1
         ENDIF
      ENDIF
 3200 CONTINUE
      IF (DEBUG.GT.0) THEN
         ASSIGN 3300 TO I99929
         GOTO 3500
      ENDIF
 3300 CONTINUE
      IF (DEBUG.GT.0) THEN
         ASSIGN 3400 TO I99975
         GOTO 1400
      ENDIF
 3400 CONTINUE
      GOTO I99977
 3500 CONTINUE
      DO 3600 I = 1, NPHI1
         DO 3550 J = 1, 2
            COMEGA(1) = COS(PHI1(J,I))
            SOMEGA(1) = SIN(PHI1(J,I))
            XX = T(1)*CTHETA(1) + T(2)*STHETA(1)*COMEGA(1) + T(3)
     +           *STHETA(1)*SOMEGA(1)
            W = (RQSUM-2*RHO(1)*XX)/(2*SIGMA(1))
            WRITE (BUFFER,9005) J, I, RR2 - XX*XX - W*W
            CALL CPRINT(BUFFER)
            WRITE (BUFFER,9006) J, I,
     +                          (RHO(2)*CTHETA(3)-(XX-RHO(1))*CTHETA(2)
     +                          -(W-SIGMA(1))*STHETA(2))
     +                          /(SIGMA(2)*STHETA(3))
            CALL CPRINT(BUFFER)
 3550    CONTINUE
 3600 CONTINUE
      GOTO I99929
 3700 CONTINUE
      IF (DEBUG.GT.0.OR.VARDBG.GT.0) THEN
         WRITE (BUFFER,9007) APHI, BPHI
         CALL CPRINT(BUFFER)
      ENDIF
      IF (BPHI.GT.APHI) THEN
         NPHI1 = 0
      ELSEIF (BPHI.GT.1.0) THEN
         NPHI1 = 0
      ELSEIF (APHI.LT.-1.0) THEN
         NPHI1 = 0
      ELSEIF (APHI.EQ.BPHI) THEN
         NPHI1 = 1
         PHI1(1,1) = ASIN(APHI)
         PHI1(2,1) = PHI1(1,1)
      ELSEIF (APHI.GT.1.0.AND.BPHI.LT.-1.0) THEN
         NPHI1 = 1
         PHI1(1,1) = -PI
         PHI1(2,1) = PI
      ELSEIF (APHI.GT.1.0.AND.BPHI.LE.1.0) THEN
         NPHI1 = 1
         PHI1(1,1) = ASIN(BPHI)
         PHI1(2,1) = PI - PHI1(1,1)
      ELSEIF (BPHI.LT.-1.0.AND.APHI.GE.-1.0) THEN
         NPHI1 = 1
         PHI1(2,1) = ASIN(APHI)
         PHI1(1,1) = -PI - PHI1(2,1)
      ELSEIF (APHI.GT.1.0.OR.BPHI.LT.-1.0) THEN
         WRITE (BUFFER,9008)
         CALL CPRINT(BUFFER)
         CALL DIE
      ELSE
         NPHI1 = 2
         IF (BPHI.LT.0.0) THEN
            PHI1(1,2) = ASIN(BPHI)
            PHI1(2,2) = ASIN(APHI)
            PHI1(1,1) = -PI - PHI1(2,2)
            PHI1(2,1) = -PI - PHI1(1,2)
         ELSE
            PHI1(1,1) = ASIN(BPHI)
            PHI1(2,1) = ASIN(APHI)
            PHI1(1,2) = PI - PHI1(2,1)
            PHI1(2,2) = PI - PHI1(1,1)
         ENDIF
      ENDIF
      DO 3800 NP = 1, NPHI1
         IF (PHI1(1,NP).GT.PHI1(2,NP)) THEN
            WRITE (BUFFER,9009) NP, PHI1(1,NP), NP, PHI1(2,NP)
            CALL CPRINT(BUFFER)
         ENDIF
 3800 CONTINUE
      GOTO I99913
 3900 CONTINUE
      RF1PHI1 = PHI1(ILIM,NP)
      KOUNT = 0
 4000 CONTINUE
      COMEGA(1) = COS(RF1PHI1)
      SOMEGA(1) = SIN(RF1PHI1)
      XX = T(1)*CTHETA(1) + T(2)*STHETA(1)*COMEGA(1) + T(3)*STHETA(1)
     +     *SOMEGA(1)
      W = (RQSUM-2*RHO(1)*XX)/(2*SIGMA(1))
      DXX = STHETA(1)*(T(3)*COMEGA(1)-T(2)*SOMEGA(1))
      DW = -RHO(1)/SIGMA(1)*DXX
      FX = (RR2-XX*XX-W*W) + 1.0E-5
      DFX = -2*(XX*DXX+W*DW)
      IF (DFX.NE.0.0) THEN
         FXDFX = FX/DFX
      ELSE
         FXDFX = 0.0
      ENDIF
      NRF1PY1 = RF1PHI1 - FXDFX
      KOUNT = KOUNT + 1
      TOOBIG = ABS(FXDFX).GT.1.0E-2 .OR. ABS(FX).GT.1.0E-3
      DONENW = KOUNT.GE.MAXNW .OR. NRF1PY1 - RF1PHI1.EQ.0.0 .OR.
     +         ABS(FXDFX).LE.EPS1 .OR. TOOBIG
      IF (DEBUG.GT.0) THEN
         WRITE (BUFFER,9010) KOUNT, NRF1PY1, RF1PHI1, FX, DFX
         CALL CPRINT(BUFFER)
      ENDIF
      IF (.NOT.(TOOBIG)) RF1PHI1 = NRF1PY1
      IF (.NOT.(DONENW)) GOTO 4000
      RF2PHI1 = PHI1(ILIM,NP)
      KOUNT = 0
 4100 CONTINUE
      COMEGA(1) = COS(RF2PHI1)
      SOMEGA(1) = SIN(RF2PHI1)
      XX = T(1)*CTHETA(1) + T(2)*STHETA(1)*COMEGA(1) + T(3)*STHETA(1)
     +     *SOMEGA(1)
      W = (RQSUM-2*RHO(1)*XX)/(2*SIGMA(1))
      DXX = STHETA(1)*(T(3)*COMEGA(1)-T(2)*SOMEGA(1))
      DW = -RHO(1)/SIGMA(1)*DXX
      FX = (RHO(2)*CTHETA(3)-(XX-RHO(1))*CTHETA(2)-(W-SIGMA(1))
     +     *STHETA(2))/(SIGMA(2)*STHETA(3))
      DFX = (-CTHETA(2)*DXX-STHETA(2)*DW)/(SIGMA(2)*STHETA(3))
      IF (FX.LT.0.0) THEN
         DFX = -DFX
         FX = -FX
      ENDIF
      FX = (FX-1.0) - 1.0E-5
      IF (DFX.NE.0.0) THEN
         FXDFX = FX/DFX
      ELSE
         FXDFX = 0.0
      ENDIF
      NRF2PY1 = RF2PHI1 - FXDFX
      KOUNT = KOUNT + 1
      TOOBIG = ABS(FXDFX).GT.1.0E-2 .OR. ABS(FX).GT.1.0E-3
      DONENW = KOUNT.GE.MAXNW .OR. NRF2PY1 - RF2PHI1.EQ.0.0 .OR.
     +         ABS(FXDFX).LE.EPS1 .OR. TOOBIG
      IF (DEBUG.GT.0) THEN
         WRITE (BUFFER,9011) KOUNT, NRF2PY1, RF2PHI1, FX, DFX
         CALL CPRINT(BUFFER)
      ENDIF
      IF (.NOT.(TOOBIG)) RF2PHI1 = NRF2PY1
      IF (.NOT.(DONENW)) GOTO 4100
      IF (RF1PHI1.GT.RF2PHI1) THEN
         TREF = RF1PHI1
         RF1PHI1 = RF2PHI1
         RF2PHI1 = TREF
      ENDIF
      GOTO I99909
 4200 CONTINUE
      DONE = .FALSE.
      FOUND = .FALSE.
      GOTO 4400
 4300 CONTINUE
      IF (DONE) THEN
         GOTO I99996
         GOTO 6100
      ENDIF
 4400 CONTINUE
      IF (NP.GT.NPHI1) GOTO 4600
      ASSIGN 4500 TO I99862
      GOTO 9900
 4500 CONTINUE
      IF (KOUNT.GT.10*BSMXIT) THEN
         WRITE (BUFFER,9012) KOUNT
         CALL CPRINT(BUFFER)
      ENDIF
      FOUND = ABS(FRTDEF).LE.EPS2 .AND. ABS(FRT).LE.EPS2
 4600 CONTINUE
      IF (FOUND) THEN
         IF (DEBUG.GT.0) THEN
            ITMAX = MAX(ITMAX,KOUNT)
            IRSTMX = MAX(IREST,IRSTMX)
            WRITE (BUFFER,9013) 'KOUNT=', KOUNT, ' ITMAX=', ITMAX,
     +                          ' IREST =', IREST, ' IRSTMX =',
     +                          IRSTMX
            CALL CPRINT(BUFFER)
         ENDIF
         DONE = .TRUE.
         ITER = ITER + 1
         ITERCT = ITERCT + 1
         ITRPHI = ITRPHI + 1
         ASSIGN 4300 TO I99856
         GOTO 15600
      ELSEIF (SGNS(1).EQ.1) THEN
         SGNS(1) = -1
         NROOT = 0
      ELSEIF (SGNS(2).EQ.1) THEN
         SGNS(2) = -1
         SGNS(1) = 1
         NROOT = 0
      ELSEIF (NP.LT.NPHI1) THEN
         ASSIGN 4700 TO I99851
         GOTO 6100
      ELSEIF (MAXDT.GT.0.AND..NOT.ADTHET.AND.ITRPHI.EQ.0) THEN
         ASSIGN 4800 TO I99848
         GOTO 6200
      ELSEIF (POMEGA(1).EQ.0.0) THEN
         ASSIGN 4900 TO I99851
         GOTO 6100
      ELSEIF (POMEGA(2).EQ.0.0) THEN
         ASSIGN 5300 TO I99851
         GOTO 6100
      ELSEIF (POMEGA(3).NE.0.0) THEN
         DONE = .TRUE.
         ITER = 0
      ELSE
         ASSIGN 5700 TO I99851
         GOTO 6100
      ENDIF
      GOTO 4300
 4700 CONTINUE
      NP = NP + 1
      SGNS(1) = 1
      SGNS(2) = 1
      NROOT = 0
      GOTO 4300
 4800 CONTINUE
      ADTHET = .TRUE.
      IF (RETRY_CLSCHN) THEN
         NP = 1
         SGNS(1) = 1
         SGNS(2) = 1
         NROOT = 0
      ENDIF
      GOTO 4300
 4900 CONTINUE
      POMEGA(1) = PI
      ADTHET = .FALSE.
      ASSIGN 5000 TO I99989
      GOTO 800
 5000 CONTINUE
      ASSIGN 5100 TO I99987
      GOTO 1000
 5100 CONTINUE
      SMLGEE = BIGGEE
      ASSIGN 5200 TO I99985
      GOTO 1600
 5200 CONTINUE
      ITRPHI = 0
      NP = 1
      SGNS(1) = 1
      SGNS(2) = 1
      NROOT = 0
      GOTO 4300
 5300 CONTINUE
      POMEGA(2) = PI
      IF (ENTZUS(1)) POMEGA(1) = 0.0
      ADTHET = .FALSE.
      ASSIGN 5400 TO I99989
      GOTO 800
 5400 CONTINUE
      ASSIGN 5500 TO I99987
      GOTO 1000
 5500 CONTINUE
      SMLGEE = BIGGEE
      ASSIGN 5600 TO I99985
      GOTO 1600
 5600 CONTINUE
      ITRPHI = 0
      NP = 1
      SGNS(1) = 1
      SGNS(2) = 1
      NROOT = 0
      GOTO 4300
 5700 CONTINUE
      POMEGA(3) = PI
      IF (ENTZUS(2)) POMEGA(2) = 0.0
      IF (ENTZUS(1)) POMEGA(1) = 0.0
      ADTHET = .FALSE.
      ASSIGN 5800 TO I99989
      GOTO 800
 5800 CONTINUE
      ASSIGN 5900 TO I99987
      GOTO 1000
 5900 CONTINUE
      SMLGEE = BIGGEE
      ASSIGN 6000 TO I99985
      GOTO 1600
 6000 CONTINUE
      ITRPHI = 0
      NP = 1
      SGNS(1) = 1
      SGNS(2) = 1
      NROOT = 0
      GOTO 4300
 
COM  OMLUPD BNJ ditch the variable format doobry for now...
COM
      IF(.NOT.(MOD(ITERCT,2) .NE. 0)) GO TO 6100
      WRITE (BUFFER,217)
  217 FORMAT('0Warning from CLSCHN -- An odd number of solutions ',
     2       'were found.')

COM      WRITE (BUFFER,217) ITERCT,((((GPHI1(ILIM,NP,IS1,IS2),ILIM=1,2),
COM     2                       IS1=1,-1,-2),IS2=1,-1,-2),NP=1,NPHI1)
COM  217 FORMAT('0Warning from CLSCHN -- An odd number of solutions ',
COM     2       'were found. ITERCT = ',I3/' GPHI1 values :'
COM     3       <4*NPHI1>(/2(1X,1PG14.7)))


      CALL CPRINT(BUFFER)
 
 6100 CONTINUE
      ITERCT = 0
      GOTO I99851
 6200 CONTINUE
      GBNDOK = .TRUE.
      IF (ABS(SMLGEE).GT.MAXG) GOTO 6800
      MODTRI = 0
      DT = MAX(0.001,MAXDT/10.0)
      ABORT = .FALSE.
      IF (ABS(SMLGEE).LE.EPS2) THEN
         WRITE (BUFFER,'(A)')
     +                      'Warning from CLSCHN - Solution was missed.'
         CALL CPRINT(BUFFER)
      ENDIF
      GOTO 6400
 6300 CONTINUE
      IF (MODDUN) GOTO 6800
 
COM
COM OMUPD DW 19/11/91 CONVERT TO MORE READABLE IF () THEN CONSTRUCT
COM
COM 99826 MODDUN = ABS(SMLGEE) .LE. EPS2 .OR.
COM     2              MODTRI .GE. MXMODT .OR.
COM     3              ABORT
COM
 
 6400 CONTINUE
      IF ((ABS(SMLGEE).LE.EPS2).OR.(MODTRI.GE.MXMODT).OR.ABORT) THEN
         MODDUN = .TRUE.
      ELSE
         MODDUN = .FALSE.
      ENDIF
 
COM
 
      IF (.NOT.(MODDUN)) THEN
         MODTRI = MODTRI + 1
         IF (NPHI1.NE.0) THEN
            SG1 = SMLGEE
            NPOLH1 = NPHI1
            DO 6420 NP = 1, NPHI1
               OLDP1(1,NP) = PHI1(1,NP)
               OLDP1(2,NP) = PHI1(2,NP)
 6420       CONTINUE
            NPOLD = NPLITL
            OLDOM1 = SMLOM1
            ASSIGN 6500 TO I99818
            GOTO 8400
         ELSE
            ASSIGN 6300 TO I99821
            GOTO 7000
         ENDIF
      ENDIF
      GOTO 6300
 6500 CONTINUE
      IF (ABS(SMLGEE).LE.EPS2) GOTO 6300
      ASSIGN 6600 TO I99815
      GOTO 9000
 6600 CONTINUE
      IF (ABS(SMLGEE).LE.EPS2) GOTO 6300
      ASSIGN 6700 TO I99985
      GOTO 1600
 6700 CONTINUE
      ASSIGN 6300 TO I99811
      GOTO 11800
 6800 CONTINUE
      RETRY_CLSCHN = ABS(SMLGEE).LE.EPS2
      IF (RETRY_CLSCHN.AND..NOT.GBNDOK) THEN
         ASSIGN 6900 TO I99985
         GOTO 1600
      ENDIF
 6900 CONTINUE
      GOTO I99848
 7000 CONTINUE
      IF (MAXG.LE.1000.0) THEN
         ASSIGN 7100 TO I99803
         GOTO 7700
      ELSE
         ASSIGN 7100 TO I99805
         GOTO 7200
      ENDIF
 7100 CONTINUE
      GOTO I99821
 7200 CONTINUE
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,*) 'Searching for G thoroughly'
         CALL CPRINT(BUFFER)
      ENDIF
      DEBUG = -DEBUG
      VARDBG = -VARDBG
      SG3 = BIGGEE
      DO 7300 I1 = -1, 1, 2
         NANGLE(1,1) = ANGLE(1,1) + I1*MAXDT
         DO 7250 I2 = -1, 1, 2
            NANGLE(2,1) = ANGLE(2,1) + I2*MAXDT
            DO 7240 I3 = -1, 1, 2
               NANGLE(3,1) = ANGLE(3,1) + I3*MAXDT
               DO 7220 I4 = -1, 1, 2
                  NANGLE(1,2) = ANGLE(1,2) + I4*MAXDT
                  DO 7215 I5 = -1, 1, 2
                     NANGLE(2,2) = ANGLE(2,2) + I5*MAXDT
                     DO 7214 I6 = -1, 1, 2
                        NANGLE(3,2) = ANGLE(3,2) + I6*MAXDT
                        DO 7212 I7 = -1, 1, 2
                           NANGLE(1,3) = ANGLE(1,3) + I7*MAXDT
                           DO 7210 I8 = -1, 1, 2
                              NANGLE(2,3) = ANGLE(2,3) + I8*MAXDT
                              DO 7208 I9 = -1, 1, 2
                                 NANGLE(3,3) = ANGLE(3,3) + I9*MAXDT
                                 ASSIGN 7202 TO I99987
                                 GOTO 1000
 7202                            CONTINUE
                                 GBNDOK = .FALSE.
                                 SGNS(1) = 1
                                 SGNS(2) = 1
                                 DO 7206 NP = 1, NPHI1
                                    OMEGA(1) = PHI1(1,NP)
     +                                     + (PHI1(2,NP)-PHI1(1,NP))/2.0
                                    ASSIGN 7204 TO I99969
                                    GOTO 14100
 7204                               CONTINUE
                                    IF (ABS(SG3).GT.ABS(GPHI)) THEN
                                       SG3 = GPHI
                                       SAVI1 = I1
                                       SAVI2 = I2
                                       SAVI3 = I3
                                       SAVI4 = I4
                                       SAVI5 = I5
                                       SAVI6 = I6
                                       SAVI7 = I7
                                       SAVI8 = I8
                                       SAVI9 = I9
                                    ENDIF
 7206                            CONTINUE
 7208                         CONTINUE
 7210                      CONTINUE
 7212                   CONTINUE
 7214                CONTINUE
 7215             CONTINUE
 7220          CONTINUE
 7240       CONTINUE
 7250    CONTINUE
 7300 CONTINUE
      DEBUG = -DEBUG
      VARDBG = -VARDBG
      IF (SG3.NE.BIGGEE) THEN
         IF (VARDBG.GT.0) THEN
            WRITE (BUFFER,9014) SAVI1, SAVI2, SAVI3, SAVI4, SAVI5,
     +                          SAVI6, SAVI7, SAVI8, SAVI9, SG3
            CALL CPRINT(BUFFER)
         ENDIF
         NANGLE(1,1) = ANGLE(1,1) + SAVI1*MAXDT
         NANGLE(2,1) = ANGLE(2,1) + SAVI2*MAXDT
         NANGLE(3,1) = ANGLE(3,1) + SAVI3*MAXDT
         NANGLE(1,2) = ANGLE(1,2) + SAVI4*MAXDT
         NANGLE(2,2) = ANGLE(2,2) + SAVI5*MAXDT
         NANGLE(3,2) = ANGLE(3,2) + SAVI6*MAXDT
         NANGLE(1,3) = ANGLE(1,3) + SAVI7*MAXDT
         NANGLE(2,3) = ANGLE(2,3) + SAVI8*MAXDT
         NANGLE(3,3) = ANGLE(3,3) + SAVI9*MAXDT
         ASSIGN 7400 TO I99987
         GOTO 1000
      ELSE
         ABORT = .TRUE.
         GOTO 7600
      ENDIF
 7400 CONTINUE
      ASSIGN 7500 TO I99985
      GOTO 1600
 7500 CONTINUE
      GBNDOK = .TRUE.
      ASSIGN 7600 TO I99811
      GOTO 11800
 7600 CONTINUE
      GOTO I99805
 7700 CONTINUE
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,*) 'Searching for G quickly'
         CALL CPRINT(BUFFER)
      ENDIF
      DO 7800 IT1 = 1, 3
         DO 7750 IT2 = 1, 3
            SAVANG = NANGLE(IT1,IT2)
            NANGLE(IT1,IT2) = ANGLE(IT1,IT2) - MAXDT
            ASSIGN 7720 TO I99780
            GOTO 8000
 7720       CONTINUE
            IF (NPHI1.GT.0) GOTO 7900
            NANGLE(IT1,IT2) = ANGLE(IT1,IT2) + MAXDT
            ASSIGN 7740 TO I99780
            GOTO 8000
 7740       CONTINUE
            IF (NPHI1.GT.0) GOTO 7900
 7750    CONTINUE
 7800 CONTINUE
      ABORT = .TRUE.
 7900 CONTINUE
      GOTO I99803
 8000 CONTINUE
      IF (SAVANG.EQ.NANGLE(IT1,IT2)) GOTO 8300
      ASSIGN 8100 TO I99987
      GOTO 1000
 8100 CONTINUE
      GBNDOK = .FALSE.
      IF (NPHI1.LE.0) GOTO 8300
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,9015) IT1, IT2, NANGLE(IT1,IT2)/DTORAD
         CALL CPRINT(BUFFER)
      ENDIF
      ASSIGN 8200 TO I99985
      GOTO 1600
 8200 CONTINUE
      GBNDOK = .TRUE.
      ASSIGN 8300 TO I99811
      GOTO 11800
 8300 CONTINUE
      GOTO I99780
 8400 CONTINUE
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,*) 'Computing secants.'
         CALL CPRINT(BUFFER)
      ENDIF
      DO 8600 IT1 = 1, 3
         DO 8550 IT2 = 1, 3
            ANG1 = NANGLE(IT1,IT2)
            ANG2 = ANG1 + DT
            IF (ANG2.LE.ANGLE(IT1,IT2)+MAXDT) THEN
               NANGLE(IT1,IT2) = ANG2
               ASSIGN 8420 TO I99987
               GOTO 1000
            ELSE
               OK = .FALSE.
               GOTO 8460
            ENDIF
 8420       CONTINUE
            GBNDOK = .FALSE.
            OK = NPHI1.GT.0
            IF (.NOT.(OK)) GOTO 8460
            ASSIGN 8440 TO I99765
            GOTO 8800
 8440       CONTINUE
            FND0G = ABS(SMLGEE).LT.EPS2 .OR.
     +              (SMLGEE*SG1.LT.0.0.AND.NPHI1.EQ.NPOLH1.AND.
     +              NPLITL.EQ.NPOLD)
            IF (.NOT.(FND0G)) GOTO 8460
            SMLGEE = 0.0
            IF (VARDBG.GT.0) THEN
               WRITE (BUFFER,*)
     +               'Secant computation stopped as solution was found.'
               CALL CPRINT(BUFFER)
            ENDIF
            GOTO 8700
 8460       CONTINUE
            IF (OK) GOTO 8520
            ANG2 = ANG1 - DT
            IF (ANG2.GE.ANGLE(IT1,IT2)-MAXDT) THEN
               NANGLE(IT1,IT2) = ANG2
               ASSIGN 8480 TO I99987
               GOTO 1000
            ELSE
               OK = .FALSE.
               GOTO 8520
            ENDIF
 8480       CONTINUE
            GBNDOK = .FALSE.
            OK = NPHI1.GT.0
            IF (.NOT.(OK)) GOTO 8520
            ASSIGN 8500 TO I99765
            GOTO 8800
 8500       CONTINUE
            FND0G = ABS(SMLGEE).LT.EPS2 .OR.
     +              (SMLGEE*SG1.LT.0.0.AND.NPHI1.EQ.NPOLH1.AND.
     +              NPLITL.EQ.NPOLD)
            IF (.NOT.(FND0G)) GOTO 8520
            SMLGEE = 0.0
            IF (VARDBG.GT.0) THEN
               WRITE (BUFFER,*)
     +               'Secant computation stopped as solution was found.'
               CALL CPRINT(BUFFER)
            ENDIF
            GOTO 8700
 8520       CONTINUE
            NANGLE(IT1,IT2) = ANG1
            IF (.NOT.(OK)) THEN
               DGDANG(IT1,IT2) = 0.0
            ELSE
               DGDANG(IT1,IT2) = (ABS(SMLGEE)-ABS(SG1))/(ANG2-ANG1)
            ENDIF
 8550    CONTINUE
 8600 CONTINUE
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,9016) ((DGDANG(IT1,IT2),IT1=1,3),IT2=1,3)
         CALL CPRINT(BUFFER)
      ENDIF
 8700 CONTINUE
      GOTO I99818
 8800 CONTINUE
      SGNS(1) = SGNSML(1)
      SGNS(2) = SGNSML(2)
      IF (NPOLH1.NE.2.OR.NPHI1.NE.2) THEN
         DENOM = OLDP1(2,NPOLH1) - OLDP1(1,1)
         IF (DENOM.NE.0.0) THEN
            P1FRAC = (OLDOM1-OLDP1(1,1))/DENOM
         ELSE
            P1FRAC = 0.0
         ENDIF
         PHIRNG = PHI1(2,NPHI1) - PHI1(1,1)
         IF (NPHI1.NE.1) THEN
            FRAC1 = (PHI1(2,1)-PHI1(1,1))/PHIRNG
            FRAC2 = (PHI1(1,2)-PHI1(1,1))/PHIRNG
            DIFF1 = P1FRAC - FRAC1
            DIFF2 = FRAC2 - P1FRAC
            IF (DIFF1*DIFF2.LE.0.0) THEN
               OMEGA(1) = PHI1(1,1) + PHIRNG*P1FRAC
            ELSEIF (DIFF1.GE.DIFF2) THEN
               OMEGA(1) = PHI1(1,2)
            ELSE
               OMEGA(1) = PHI1(2,1)
            ENDIF
            IF (OMEGA(1).GT.PHI1(2,1)) THEN
               NPNEW = 2
            ELSE
               NPNEW = 1
            ENDIF
         ELSE
            OMEGA(1) = PHI1(1,1) + PHIRNG*P1FRAC
            NPNEW = 1
         ENDIF
      ELSE
         P1FRAC = (OLDOM1-OLDP1(1,NPOLD))
     +            /(OLDP1(2,NPOLD)-OLDP1(1,NPOLD))
         PHIRNG = PHI1(2,NPOLD) - PHI1(1,NPOLD)
         OMEGA(1) = PHI1(1,NPOLD) + P1FRAC*PHIRNG
         NPNEW = NPOLD
      ENDIF
      ASSIGN 8900 TO I99969
      GOTO 14100
 8900 CONTINUE
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,9017) OMEGA(1), GPHI
         CALL CPRINT(BUFFER)
      ENDIF
      SMLGEE = GPHI
      NPLITL = NPNEW
      GOTO I99765
 9000 CONTINUE
      DEN = 0.0
      DO 9100 IT1 = 1, 3
         DO 9050 IT2 = 1, 3
            DEN = DEN + DGDANG(IT1,IT2)**2
            OLDANG(IT1,IT2) = NANGLE(IT1,IT2)
 9050    CONTINUE
 9100 CONTINUE
      IF (DEN.GE.EPS1) THEN
         SMLGEE = SG1
         NPOLH1 = NPHI1
         NPOLD = NPLITL
         LENSD = SQRT(DEN)
         STEPSD = 1.0
         NSDSTP = 0
         GOTO 9300
      ELSE
         ABORT = .TRUE.
         GOTO I99815
         GOTO 9900
      ENDIF
 9200 CONTINUE
      IF (STOPSD) THEN
         ABORT = LENMV.LT.EPS1
         GOTO I99815
         GOTO 9900
      ENDIF
 9300 CONTINUE
      NSDSTP = NSDSTP + 1
      OLDG = SMLGEE
      SGNREF(1) = SGNSML(1)
      SGNREF(2) = SGNSML(2)
      REFOM1 = SMLOM1
      REFNP = NPLITL
      LNSHFT = 0.0
      LENMV = 0.0
      DO 9400 IT1 = 1, 3
         DO 9350 IT2 = 1, 3
            PRVANG(IT1,IT2) = NANGLE(IT1,IT2)
            NANGLE(IT1,IT2) = OLDANG(IT1,IT2) - STEPSD*DGDANG(IT1,IT2)
     +                        *ABS(SG1)/DEN
            IF (NANGLE(IT1,IT2).LT.ANGLE(IT1,IT2)-MAXDT) THEN
               NANGLE(IT1,IT2) = ANGLE(IT1,IT2) - MAXDT
            ELSEIF (NANGLE(IT1,IT2).GT.ANGLE(IT1,IT2)+MAXDT) THEN
               NANGLE(IT1,IT2) = ANGLE(IT1,IT2) + MAXDT
            ENDIF
            NMSHFT(IT1,IT2) = NANGLE(IT1,IT2) - PRVANG(IT1,IT2)
            LNSHFT = LNSHFT + NMSHFT(IT1,IT2)**2
            LENMV = LENMV + (NANGLE(IT1,IT2)-OLDANG(IT1,IT2))**2
            IF (VARDBG.GT.0) THEN
               WRITE (BUFFER,9018) IT1, IT2, NMSHFT(IT1,IT2)
               CALL CPRINT(BUFFER)
            ENDIF
 9350    CONTINUE
 9400 CONTINUE
      ASSIGN 9500 TO I99987
      GOTO 1000
 9500 CONTINUE
      GBNDOK = .FALSE.
      IF (NPHI1.NE.0) THEN
         ASSIGN 9600 TO I99765
         GOTO 8800
      ELSE
         SMLGEE = BIGGEE
      ENDIF
 9600 CONTINUE
      FND0G = ABS(SMLGEE).LT.EPS2 .OR.
     +        (SMLGEE*SG1.LT.0.0.AND.NPHI1.EQ.NPOLH1.AND.
     +        NPLITL.EQ.NPOLD)
      IF (.NOT.(FND0G)) THEN
         LENMV = SQRT(LENMV)
         LNSHFT = SQRT(LNSHFT)
         STEPSD = STEPSD*1.5
         IF (VARDBG.GT.0) THEN
            WRITE (BUFFER,'(A,1PG12.5,0PF8.5)') ' LNSHFT,LENMV =',
     +             LNSHFT, LENMV
            CALL CPRINT(BUFFER)
         ENDIF
         IF (ABS(SMLGEE).LE.ABS(OLDG).OR.NSDSTP.LE.1) GOTO 9800
         DO 9650 IT1 = 1, 3
            DO 9620 IT2 = 1, 3
               NANGLE(IT1,IT2) = PRVANG(IT1,IT2)
 9620       CONTINUE
 9650    CONTINUE
         ASSIGN 9700 TO I99987
         GOTO 1000
      ELSE
         SMLGEE = 0.0
         IF (VARDBG.GT.0) THEN
            WRITE (BUFFER,*) 'SD stopping as solution found.'
            CALL CPRINT(BUFFER)
         ENDIF
         GOTO 9800
      ENDIF
 9700 CONTINUE
      GBNDOK = .FALSE.
      SMLGEE = OLDG
      SGNSML(1) = SGNREF(1)
      SGNSML(2) = SGNREF(2)
      SMLOM1 = REFOM1
      NPLITL = REFNP
 9800 CONTINUE
      DRPSTP = (ABS(SMLGEE).GE.ABS(OLDG).OR.NPHI1.EQ.0) .AND.
     +         NSDSTP.EQ.1 .AND. LNSHFT.GE.EPS1
      IF (.NOT.(DRPSTP)) THEN
         STOPSD = ABS(SMLGEE).LT.EPS2 .OR. ABS(SMLGEE).GE.ABS(OLDG) .OR.
     +            LNSHFT.LT.EPS1
      ELSE
         STEPSD = 0.1
         STOPSD = .FALSE.
      ENDIF
      GOTO 9200
 9900 CONTINUE
      MULLER = .TRUE.
      KOUNT = 0
      IREST = 0
      LOWPHI = PHI1(1,NP)
      HIPHI = PHI1(2,NP)
      DIFPHI = HIPHI - LOWPHI
      GLOPHI = GPHI1(1,NP,SGNS(1),SGNS(2))
      GHIPHI = GPHI1(2,NP,SGNS(1),SGNS(2))
      IF (NROOT.EQ.0) THEN
         IF (ABS(GLOPHI).LE.EPS2) THEN
            RT = -10.0
            NROOT = 1
            OMEGA(1) = LOWPHI
            FRTDEF = GLOPHI
            FRT = GLOPHI
            IF (.TRUE.) GOTO 11200
         ELSEIF (ABS(GHIPHI).LE.EPS2) THEN
            RT = 10.0
            NROOT = 1
            OMEGA(1) = HIPHI
            FRTDEF = GHIPHI
            FRT = GHIPHI
            IF (.TRUE.) GOTO 11200
         ENDIF
      ENDIF
      NROOT = NROOT + 1
      PRODG = GLOPHI*GHIPHI
      IF (PRODG.GE.0.0.OR.NROOT.NE.1) THEN
         IF (DEBUG.GT.0) THEN
            WRITE (BUFFER,9019)
            CALL CPRINT(BUFFER)
         ENDIF
         POSSEN = .FALSE.
         NEGSEN = .FALSE.
         RT = -10.0
         ASSIGN 10000 TO I99703
         GOTO 13400
      ELSE
         IF (DEBUG.GT.0) THEN
            WRITE (BUFFER,9020)
            CALL CPRINT(BUFFER)
         ENDIF
         FRF = GLOPHI
         GRF = GHIPHI
         ARF = LOWPHI
         BRF = HIPHI
         ASSIGN 11300 TO I99705
         GOTO 11500
      ENDIF
      GOTO 11300
10000 CONTINUE
      RT = 10.0
      ASSIGN 10100 TO I99703
      GOTO 13400
10100 CONTINUE
      MAXIT = BSMXIT
      RTS(NROOT) = 0.0
10200 CONTINUE
      H = 0.5
      RT = RTS(NROOT) + H
      ASSIGN 10300 TO I99703
      GOTO 13400
10300 CONTINUE
      ASSIGN 10400 TO I99699
      GOTO 11400
10400 CONTINUE
      DELFPR = FRTDEF
      RT = RTS(NROOT) - H
      ASSIGN 10500 TO I99703
      GOTO 13400
10500 CONTINUE
      ASSIGN 10600 TO I99699
      GOTO 11400
10600 CONTINUE
      FRTPRV = FRTDEF
      DELFPR = FRTPRV - DELFPR
      RT = RTS(NROOT)
      ASSIGN 10700 TO I99703
      GOTO 13400
10700 CONTINUE
      ASSIGN 10800 TO I99699
      GOTO 11400
10800 CONTINUE
      LAMBDA = -0.5
10900 CONTINUE
      DELF = FRTDEF - FRTPRV
      DFPRLM = DELFPR*LAMBDA
      NUM = -FRTDEF*(1.0+LAMBDA)*2.0
      G = (1.0+LAMBDA*2.0)*DELF - LAMBDA*DFPRLM
      SQR = G*G + 2.0*NUM*LAMBDA*(DELF-DFPRLM)
      IF (SQR.LT.0.0) SQR = 0.0
      SQR = SQRT(SQR)
      DEN = G + SQR
      IF (G*SQR.LT.0) DEN = G - SQR
      IF (ABS(DEN).EQ.0) DEN = 1.0
      LAMBDA = NUM/DEN
      FRTPRV = FRTDEF
      DELFPR = DELF
      IF (ABS(LAMBDA).GT.100) LAMBDA = SIGN(100.0,LAMBDA)
      H = H*LAMBDA
      RT = RT + H
      IF (KOUNT.LE.MAXIT) THEN
         ASSIGN 11000 TO I99703
         GOTO 13400
      ELSEIF (KOUNT.LE.10*BSMXIT) THEN
         CARYON = (POSSEN.AND.NEGSEN) .OR.
     +            (MOD(NROOT,2).EQ.0.AND.PRODG.GT.0.0) .OR.
     +            (MOD(NROOT,2).EQ.1.AND.PRODG.LT.0.0)
         IF (POSSEN.AND.NEGSEN) THEN
            ARF = OMEGA1POS
            BRF = OMEGA1NEG
            FRF = FRTDEFPOS
            GRF = FRTDEFNEG
            ASSIGN 11200 TO I99705
            GOTO 11500
         ELSEIF (CARYON) THEN
            IF (DEBUG.GT.0) THEN
               WRITE (BUFFER,9021) MAXIT
               CALL CPRINT(BUFFER)
            ENDIF
            MAXIT = MAXIT + BSMXIT
            ASSIGN 11000 TO I99703
            GOTO 13400
         ENDIF
      ENDIF
      GOTO 11200
11000 CONTINUE
      IF (AMAX1(ABS(FRT),ABS(FRTDEF)).LE.EPS2) GOTO 11200
      ASSIGN 11100 TO I99699
      GOTO 11400
11100 CONTINUE
      OVERLP = AMAX1(ABS(OMEGA(1)-OMPRV1),ABS(OMPRV1-OMPRV2)).EQ.0.0
      IF (.NOT.(OVERLP)) THEN
         IF (ABS(FRTDEF).LE.10.0*ABS(FRTPRV)) GOTO 10900
         H = H/2.0
         LAMBDA = LAMBDA/2.0
         RT = RT - H
         IF (.TRUE.) THEN
            ASSIGN 11000 TO I99703
            GOTO 13400
         ENDIF
      ELSE
         IF (RTS(NROOT).LE.0.0) THEN
            RTS(NROOT) = 0.6 - RTS(NROOT)
         ELSE
            RTS(NROOT) = -RTS(NROOT)
         ENDIF
         IREST = IREST + 1
         GOTO 10200
      ENDIF
11200 CONTINUE
      RTS(NROOT) = RT
      PHIRTS(NROOT) = OMEGA(1)
      IF (DEBUG.GT.0) THEN
         WRITE (BUFFER,9022) OMEGA(1), FRT, KOUNT
         CALL CPRINT(BUFFER)
      ENDIF
11300 CONTINUE
      GOTO I99862
11400 CONTINUE
      IF (FRTDEFPREV*FRTDEF.GE.0.0.OR.KOUNT.LE.BSMXIT) THEN
         GOTO I99699
      ELSE
         IF (DEBUG.GT.0) THEN
            WRITE (BUFFER,*) 'Switching to modified regula falsi.'
            CALL CPRINT(BUFFER)
         ENDIF
         FRF = FRTDEF
         GRF = FRTDEFPREV
         ARF = OMEGA(1)
         BRF = OMPRV1
         ASSIGN 11200 TO I99705
      ENDIF
11500 CONTINUE
      MULLER = .FALSE.
      MAXIT = KOUNT + MXITRF
      WRF = ARF
      FRTDEF = FRF
      SFRF = SIGN(1.0,FRF)
11600 CONTINUE
      WRF = (FRF*BRF-GRF*ARF)/(FRF-GRF)
      SWRF = SIGN(1.0,FRTDEF)
      OMEGA(1) = WRF
      RT = WRF
      ASSIGN 11700 TO I99678
      GOTO 13600
11700 CONTINUE
      IF (SFRF*FRTDEF.GE.0.0) THEN
         ARF = WRF
         FRF = FRTDEF
         IF (SWRF*FRTDEF.GT.0.0) GRF = GRF/2.0
      ELSE
         BRF = WRF
         GRF = FRTDEF
         IF (SWRF*FRTDEF.GT.0.0) FRF = FRF/2.0
      ENDIF
      DONERF = KOUNT.GT.MAXIT .OR. MAX(ABS(FRT),ABS(FRTDEF))
     +         .LE.EPS2 .OR. ABS(BRF-ARF).LE.EPS1
      IF (.NOT.(DONERF)) GOTO 11600
      PHIRTS(NROOT) = RT
      IF (ABS(FRT).GE.EPS1) THEN
         RATIO = FRTDEF/FRT
      ELSE
         RATIO = FRTDEF/EPS1
      ENDIF
      IF (ABS(BRF-ARF).LE.EPS1.AND.ABS(RATIO).LT.100.0) THEN
         FRTDEF = EPS2
         FRT = EPS2
      ENDIF
C
C
COM OMUPD BNJ update code for erfi function...
C
      RTEMP1 = 2.0*(RT-LOWPHI)/(HIPHI-LOWPHI)-1.0
      CALL ABMERFI(RTEMP1,RTS(NROOT),IER)
C
      IF (DEBUG.GT.1) THEN
         WRITE (BUFFER,9023) RTS(NROOT), IER
         CALL CPRINT(BUFFER)
      ENDIF
      GOTO I99705
11800 CONTINUE
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,*) 'Finding minimum G'
         CALL CPRINT(BUFFER)
      ENDIF
      SMLGEE = BIGGEE
      DO 11900 NP = 1, NPHI1
         DO 11850 IS2 = 1, -1, -2
            SGNS(2) = IS2
            DO 11820 IS1 = 1, -1, -2
               SGNS(1) = IS1
               ASSIGN 11820 TO I99668
               GOTO 12000
11820       CONTINUE
11850    CONTINUE
11900 CONTINUE
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,*) 'Smallest G = ', SMLGEE
         CALL CPRINT(BUFFER)
      ENDIF
      GOTO I99811
12000 CONTINUE
      MULLER = .TRUE.
      KOUNT = 0
      LOWPHI = PHI1(1,NP)
      HIPHI = PHI1(2,NP)
      DIFPHI = HIPHI - LOWPHI
      GLOPHI = GPHI1(1,NP,SGNS(1),SGNS(2))
      GHIPHI = GPHI1(2,NP,SGNS(1),SGNS(2))
      GPHI = GLOPHI
      OMEGA(1) = LOWPHI
      ASSIGN 12100 TO I99666
      GOTO 14800
12100 CONTINUE
      GPHI = GHIPHI
      OMEGA(1) = HIPHI
      ASSIGN 12200 TO I99666
      GOTO 14800
12200 CONTINUE
      IF (ABS(SMLGEE).LT.EPS2) GOTO 13200
      NROOT = 1
      PRODG = GLOPHI*GHIPHI
      IF (PRODG.GE.0.0) THEN
         POSSEN = .FALSE.
         NEGSEN = .FALSE.
         MAXIT = MXITMN
         RTS(NROOT) = 0.0
         H = 0.5
         RT = RTS(NROOT) + H
         ASSIGN 12300 TO I99662
         GOTO 13900
      ELSE
         IF (VARDBG.GT.0) THEN
            WRITE (BUFFER,9024)
            CALL CPRINT(BUFFER)
         ENDIF
         SMLGEE = 0.0
         GOTO 13200
      ENDIF
12300 CONTINUE
      ASSIGN 12400 TO I99660
      GOTO 13300
12400 CONTINUE
      DELFPR = GPHI
      RT = RTS(NROOT) - H
      ASSIGN 12500 TO I99662
      GOTO 13900
12500 CONTINUE
      ASSIGN 12600 TO I99660
      GOTO 13300
12600 CONTINUE
      FRTPRV = GPHI
      DELFPR = FRTPRV - DELFPR
      RT = RTS(NROOT)
      ASSIGN 12700 TO I99662
      GOTO 13900
12700 CONTINUE
      ASSIGN 12800 TO I99660
      GOTO 13300
12800 CONTINUE
      LAMBDA = -0.5
12900 CONTINUE
      DELF = GPHI - FRTPRV
      DFPRLM = DELFPR*LAMBDA
      NUM = -GPHI*(1.0+LAMBDA)*2.0
      G = (1.0+LAMBDA*2.0)*DELF - LAMBDA*DFPRLM
      SQR = G*G + 2.0*NUM*LAMBDA*(DELF-DFPRLM)
      IF (SQR.LT.0.0) SQR = 0.0
      SQR = SQRT(SQR)
      DEN = G + SQR
      IF (G*SQR.LT.0) DEN = G - SQR
      IF (ABS(DEN).EQ.0) DEN = 1.0
      LAMBDA = NUM/DEN
      FRTPRV = GPHI
      DELFPR = DELF
      IF (ABS(LAMBDA).GT.100) LAMBDA = SIGN(100.0,LAMBDA)
      H = H*LAMBDA
      RT = RT + H
      IF (KOUNT.GT.MAXIT) GOTO 13200
      ASSIGN 13000 TO I99662
      GOTO 13900
13000 CONTINUE
      IF (ABS(SMLGEE).LE.EPS2) GOTO 13200
      ASSIGN 13100 TO I99660
      GOTO 13300
13100 CONTINUE
      OVERLP = AMAX1(ABS(OMEGA(1)-OMPRV1),ABS(OMPRV1-OMPRV2)).EQ.0.0
      IF (.NOT.(OVERLP)) THEN
         IF (ABS(GPHI).LE.10.0*ABS(FRTPRV)) GOTO 12900
         H = H/2.0
         LAMBDA = LAMBDA/2.0
         RT = RT - H
         IF (.TRUE.) THEN
            ASSIGN 13000 TO I99662
            GOTO 13900
         ENDIF
      ENDIF
13200 CONTINUE
      IF (VARDBG.GT.0) THEN
         WRITE (BUFFER,9025) OMEGA(1), SMLGEE, KOUNT
         CALL CPRINT(BUFFER)
      ENDIF
      GOTO I99668
13300 CONTINUE
      IF (GPPREV*GPHI.GE.0.0) THEN
         GOTO I99660
      ELSE
         IF (VARDBG.GT.0) THEN
            WRITE (BUFFER,*) 'Sign change seen so zero inferred.'
            CALL CPRINT(BUFFER)
         ENDIF
         SMLGEE = 0
         GOTO 13200
      ENDIF
13400 CONTINUE
      OMPRV2 = OMPRV1
      OMPRV1 = OMEGA(1)
      OMEGA(1) = LOWPHI + DIFPHI*0.5*(1.0+ABMERF(RT))
      FRTDEFPREV = FRTDEF
      ASSIGN 13500 TO I99678
      GOTO 13600
13500 CONTINUE
      IF (FRTDEF.GT.0) THEN
         POSSEN = .TRUE.
         FRTDEFPOS = FRTDEF
         OMEGA1POS = OMEGA(1)
      ELSEIF (FRTDEF.LT.0) THEN
         NEGSEN = .TRUE.
         FRTDEFNEG = FRTDEF
         OMEGA1NEG = OMEGA(1)
      ENDIF
      GOTO I99703
13600 CONTINUE
      KOUNT = KOUNT + 1
      IF (OMEGA(1).EQ.LOWPHI) THEN
         GPHI = GLOPHI
      ELSEIF (OMEGA(1).NE.HIPHI) THEN
         ASSIGN 13700 TO I99969
         GOTO 14100
      ELSE
         GPHI = GHIPHI
      ENDIF
13700 CONTINUE
      FRT = GPHI
      FRTDEF = FRT
      DO 13800 J = 1, NROOT - 1
         DEN = OMEGA(1) - PHIRTS(J)
         IF (ABS(DEN).GT.EPS1) THEN
            FRTDEF = FRTDEF/DEN
         ELSEIF (.NOT.(MULLER)) THEN
            DEN = EPS1
            FRTDEF = FRTDEF/DEN
         ELSE
            IF (ABS(RT).LE.3.2) THEN
               RTS(NROOT) = RT + .001
            ELSE
               RTS(NROOT) = 2*(INT(RT)-RT)
            ENDIF
            IF (DEBUG.GT.1) THEN
               WRITE (BUFFER,9026) RT, FRT
               CALL CPRINT(BUFFER)
            ENDIF
            IF (KOUNT.LE.11.0*BSMXIT) GOTO 10200
            GOTO 11200
         ENDIF
13800 CONTINUE
      IF (DEBUG.GT.1.OR.VARDBG.GT.1) THEN
         WRITE (BUFFER,9026) RT, OMEGA(1), FRT, FRTDEF
         CALL CPRINT(BUFFER)
      ENDIF
      GOTO I99678
13900 CONTINUE
      KOUNT = KOUNT + 1
      GPPREV = GPHI
      OMEGA(1) = LOWPHI + DIFPHI*0.5*(1.0+ABMERF(RT))
      IF (OMEGA(1).EQ.LOWPHI) THEN
         GPHI = GLOPHI
      ELSEIF (OMEGA(1).NE.HIPHI) THEN
         ASSIGN 14000 TO I99969
         GOTO 14100
      ELSE
         GPHI = GHIPHI
      ENDIF
14000 CONTINUE
      GOTO I99662
14100 CONTINUE
      ASSIGN 14200 TO I99632
      GOTO 14900
14200 CONTINUE
      ASSIGN 14300 TO I99630
      GOTO 15100
14300 CONTINUE
      ASSIGN 14400 TO I99628
      GOTO 15200
14400 CONTINUE
      ASSIGN 14500 TO I99626
      GOTO 15300
14500 CONTINUE
      ASSIGN 14600 TO I99624
      GOTO 15400
14600 CONTINUE
      ASSIGN 14700 TO I99666
      GOTO 14800
14700 CONTINUE
      GOTO I99969
14800 CONTINUE
      IF (ABS(GPHI).LT.ABS(SMLGEE)) THEN
         SMLGEE = GPHI
         SMLOM1 = OMEGA(1)
         SGNSML(1) = SGNS(1)
         SGNSML(2) = SGNS(2)
         NPLITL = NP
      ENDIF
      GOTO I99666
14900 CONTINUE
      DO 15000 I = 1, 3
         TA(I) = S(I) - Q(I,0)
15000 CONTINUE
      R(1) = CTHETA(0)*TA(1) + STHETA(0)*TA(2)
      R(2) = -STHETA(0)*TA(1) + CTHETA(0)*TA(2)
      R(3) = TA(3)
      COMEGA(1) = COS(OMEGA(1))
      SOMEGA(1) = SIN(OMEGA(1))
      TA(1) = R(1)
      TA(2) = COMEGA(1)*R(2) + SOMEGA(1)*R(3)
      TA(3) = -SOMEGA(1)*R(2) + COMEGA(1)*R(3)
      R(1) = CTHETA(1)*TA(1) + STHETA(1)*TA(2)
      R(2) = -STHETA(1)*TA(1) + CTHETA(1)*TA(2)
      R(3) = TA(3)
      XX = R(1)
      YY = R(2)
      ZZ = R(3)
      GOTO I99632
15100 CONTINUE
      W = (RQSUM-2*RHO(1)*XX)/(2*SIGMA(1))
      DENOM = RR2 - XX*XX
      SQR = RR2 - XX*XX - W*W
      IF (SQR.LT.0.0) THEN
         IF (ABS(SQR).GT.1.0E-3) THEN
            WRITE (BUFFER,'(A,1PG14.7)')
     +              ' Warning from CLSCHN -- Omega2 square is ', SQR
            CALL CPRINT(BUFFER)
         ENDIF
         SQR = 0.0
      ENDIF
      SQR = SQRT(SQR)
      COMEGA(2) = (YY*W+SGNS(1)*ZZ*SQR)/DENOM
      SOMEGA(2) = (ZZ*W-SGNS(1)*YY*SQR)/DENOM
      Q3(1) = XX - RHO(1)
      Q3(2) = W - SIGMA(1)
      Q3(3) = SGNS(1)*SQR
      GOTO I99630
15200 CONTINUE
      COMEGA(4) = (RHO(2)*CTHETA(3)-Q3(1)*CTHETA(2)-Q3(2)*STHETA(2))
     +            /(SIGMA(2)*STHETA(3))
      SQR = 1.0 - COMEGA(4)**2
      IF (SQR.LT.0.0) THEN
         IF (ABS(SQR).GT.1.0E-3) THEN
            WRITE (BUFFER,'(A,1PG14.7)')
     +              ' Warning from CLSCHN -- Omega4 square is ', SQR
            CALL CPRINT(BUFFER)
         ENDIF
         SQR = 0.0
         COMEGA(4) = SIGN(1.0,COMEGA(4))
      ENDIF
      SQR = SQRT(SQR)
      SOMEGA(4) = SGNS(2)*SQR
      GOTO I99628
15300 CONTINUE
      A = RHO(2)*STHETA(3) + SIGMA(2)*COMEGA(4)*CTHETA(3)
      B = Q3(2)*CTHETA(2) - Q3(1)*STHETA(2)
      C = SIGMA(2)*SOMEGA(4)
      DENOM = B**2 + Q3(3)**2
      COMEGA(3) = (A*B+Q3(3)*C)/DENOM
      SOMEGA(3) = (A*Q3(3)-B*C)/DENOM
      GOTO I99626
15400 CONTINUE
      TA(1) = CTHETA(4)
      TA(2) = STHETA(4)
      TA(3) = 0.0
      DO 15500 I = 4, 1, -1
         TB(1) = TA(1)
         TB(2) = COMEGA(I)*TA(2) - SOMEGA(I)*TA(3)
         TB(3) = SOMEGA(I)*TA(2) + COMEGA(I)*TA(3)
         TA(1) = CTHETA(I-1)*TB(1) - STHETA(I-1)*TB(2)
         TA(2) = STHETA(I-1)*TB(1) + CTHETA(I-1)*TB(2)
         TA(3) = TB(3)
15500 CONTINUE
      GPHI = DOT(U,TA,3) - CTHETA(5)
      GOTO I99624
15600 CONTINUE
      IF (OMEGA(1).EQ.LOWPHI.OR.OMEGA(1).EQ.HIPHI) THEN
         ASSIGN 15700 TO I99969
         GOTO 14100
      ENDIF
15700 CONTINUE
      OMEGA(2) = ATAN2(SOMEGA(2),COMEGA(2))
      OMEGA(3) = ATAN2(SOMEGA(3),COMEGA(3))
      OMEGA(4) = ATAN2(SOMEGA(4),COMEGA(4))
      ASSIGN 15800 TO I99612
      GOTO 16300
15800 CONTINUE
      ASSIGN 15900 TO I99610
      GOTO 16500
15900 CONTINUE
      ASSIGN 16000 TO I99608
      GOTO 16800
16000 CONTINUE
      DO 16100 I = 1, 5, 2
         IF (POMEGA((I+1)/2).NE.0.0) THEN
            OMEGA(I) = OMEGA(I) + PI
         ENDIF
16100 CONTINUE
      DO 16200 I = 1, 6
         OMEGA(I) = OMEGA(I) - 2*PI*NINT(OMEGA(I)/(2*PI))
16200 CONTINUE
      GOTO I99856
16300 CONTINUE
      TA(1) = CTHETA(0)*U(1) + STHETA(0)*U(2)
      TA(2) = -STHETA(0)*U(1) + CTHETA(0)*U(2)
      TA(3) = U(3)
      DO 16400 I = 1, 4
         ASSIGN 16400 TO I99602
         GOTO 16700
16400 CONTINUE
      IF (DEBUG.GT.0) THEN
         WRITE (BUFFER,9027) TA(1), CTHETA(5)
         CALL CPRINT(BUFFER)
      ENDIF
      COMEGA(5) = TA(2)/STHETA(5)
      SOMEGA(5) = TA(3)/STHETA(5)
      OMEGA(5) = ATAN2(SOMEGA(5),COMEGA(5))
      GOTO I99612
16500 CONTINUE
      TA(1) = CTHETA(0)*V(1) + STHETA(0)*V(2)
      TA(2) = -STHETA(0)*V(1) + CTHETA(0)*V(2)
      TA(3) = V(3)
      DO 16600 I = 1, 5
         ASSIGN 16600 TO I99602
         GOTO 16700
16600 CONTINUE
      IF (DEBUG.GT.0) THEN
         WRITE (BUFFER,9028) TA(1)
         CALL CPRINT(BUFFER)
      ENDIF
      COMEGA(6) = TA(2)
      SOMEGA(6) = TA(3)
      OMEGA(6) = ATAN2(SOMEGA(6),COMEGA(6))
      GOTO I99610
16700 CONTINUE
      TB(1) = TA(1)
      TB(2) = COMEGA(I)*TA(2) + SOMEGA(I)*TA(3)
      TB(3) = -SOMEGA(I)*TA(2) + COMEGA(I)*TA(3)
      TA(1) = CTHETA(I)*TB(1) + STHETA(I)*TB(2)
      TA(2) = -STHETA(I)*TB(1) + CTHETA(I)*TB(2)
      TA(3) = TB(3)
      GOTO I99602
16800 CONTINUE
      DO 16900 I = 6, 0, -1
         IF (I.NE.0) THEN
            NEWX(I) = 0.0
            NEWY(I) = 0.0
            NEWZ(I) = 0.0
         ENDIF
         DO 16850 J = I + 1, 6
            TA(1) = NEWX(J)
            TA(2) = NEWY(J)*COMEGA(I+1) - NEWZ(J)*SOMEGA(I+1)
            TA(3) = NEWY(J)*SOMEGA(I+1) + NEWZ(J)*COMEGA(I+1)
            NEWX(J) = TA(1)*CTHETA(I) - TA(2)*STHETA(I) + P(1,I)
            NEWY(J) = TA(1)*STHETA(I) + TA(2)*CTHETA(I) + P(2,I)
            NEWZ(J) = TA(3) + P(3,I)
16850    CONTINUE
16900 CONTINUE
      IND = ATMIND(1,1)
      DO 17000 I = 1, 6
         TA(1) = X(IND)
         TA(2) = Y(IND)
         TA(3) = Z(IND)
         DO 16950 J = 1, 3
            TA(J) = TA(J) + NEWX(I)*U0(J) + NEWY(I)*V0(J) + NEWZ(I)
     +              *W0(J)
16950    CONTINUE
         NEWX(I) = TA(1)
         NEWY(I) = TA(2)
         NEWZ(I) = TA(3)
17000 CONTINUE
      GOTO I99608
C
C OMUPD rkw 13/08/92 swapped inappropriate '/' with ','
C
C9001 FORMAT ('0Warning from CLSCHN -- Boundaries failed to match.'/
C    +        ' ILIM = ',I1,' NP = ',I1/' G(PHI1) =',4(1X,1PG14.7))
 9001 FORMAT ('0Warning from CLSCHN -- Boundaries failed to match.',
     +        ' ILIM = ',I1,' NP = ',I1,' G(PHI1) =',4(1X,1PG14.7))
C
 9002 FORMAT (' RR,QR1,QR2 = ',1P3G14.7)
 9003 FORMAT (' LB,UB = ',1P2G14.7)
 9004 FORMAT (' PHI1: ',4(1PG14.7,1X))
 9005 FORMAT (' Bound for eq. 28 for PHI1(',I1,',',I1,') is ',1PG14.7)
 9006 FORMAT (' Bound for eq. 34 for PHI1(',I1,',',I1,') is ',1PG14.7)
 9007 FORMAT (' APHI and BPHI for GET-PHI-RANGES are ',1P2G14.7)
 9008 FORMAT ('0Error in CLSCHN -- Fell through ','PHI CONDITIONAL')
 9009 FORMAT (' Phi1 values out of order PHI1(1,',I1,') = ',1PG14.7,
     +        ' PHI1(2,',I1,') = ',1PG14.7)
 9010 FORMAT (' Eq. 28 bounds:',I2,':',1P4G14.7)
 9011 FORMAT (' Eq. 34 bounds:',I2,':',1P4G14.7)
 9012 FORMAT ('0KOUNT exceeded 10*BSMXIT.  KOUNT = ',I12)
 9013 FORMAT (1X,4(A,I3))
 9014 FORMAT (' Search for G, Angle sign = ',9I3,' G =',1PG14.7)
 9015 FORMAT (' G found for NANGLE(',I1,',',I1,')=',F7.3)
 9016 FORMAT (' DG/DANGLE =',3(/1X,1P3G14.7))
 9017 FORMAT (' G(',F7.4,') =',1PG14.7)
 9018 FORMAT (' Shift for ANGLE(',I1,',',I1,')=',1PG14.7)
 9019 FORMAT ('Muller''s method is being used to solve the eqn.')
 9020 FORMAT (' Regula falsi being used to solve the equation.')
 9021 FORMAT (' Incrementing MAXIT to ',I4,' as more roots are likely.')
 9022 FORMAT (' A Root is ',1PG14.7,' F is ',1PG14.7,I5,' iterations')
 9023 FORMAT (' A solution in unconstrained space is ',1PG14.7,
     +        '  IER = ',I3)
 9024 FORMAT (' Sign change indicates G has zero value.')
 9025 FORMAT (' A Minimum G is at ',1PG14.7,' G is ',1PG14.7,I5,
     +        ' iterations')
 9026 FORMAT (1X,1P4G14.7)
 9027 FORMAT (1X,1PG14.7,'should == COS(THETA(5)) = ',1PG14.7)
 9028 FORMAT (' The following number should be zero ',1PG14.7)
      END



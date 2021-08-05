C
C     Alternate implementation of CARTX which performs the basic BILDER
C     function, the construction of an atom's position from the three
C     atoms previous on the chain using the bond length, bond angle, and
C     torsion angle. We construct the position of atom I4 from I1, I2,
C     and I3. The bond is between I3 and I4, the angle is between I2,
C     I3, and I4, and the torsion is through all four. X, Y, Z are the
C     coordinate arrays. This implementation is faster than CARTX.
C
C     The basic algorithm is to translate and rotate the atoms so that
C     X3 is on the origin, X2 on is the negative x axis, and X1 is on
C     the XY plane. The position of X4 is easily generated, and then the
C     inverse transformation is performed to bring X4 to its proper
C     place.
C
COM OMUPD BNJ 2/9/91

      SUBROUTINE CARTX2(X,Y,Z,I1,I2,I3,I4,BOND,THETA,PHI)
C
      IMPLICIT NONE
C
C -- Includes:-
C
      include "values.inc"
C
C -- Declarations:-
C
      REAL X(*), Y(*), Z(*)
      INTEGER I1, I2, I3, I4
      REAL BOND, THETA, PHI
      REAL CTHT, STHT, BSIN, CPHI, SPHI
      REAL X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4
      REAL LYZ1, OVLYZ1, YY4, ZZ4, Y1O, Z1O, LXZ22, L2, LXZ2, OVL2,
     +     OVLXZ2
      REAL X2O, Z2O, XZ2O, Y2O, XX1, XX4
C
      CHARACTER*100 BUFFER
C
      REAL ETA , ETA2
      PARAMETER ( ETA = 1.0E-7, ETA2 = ETA*ETA )
C
C -- Code:-
C
      IF (X(I1).EQ.ANUM.OR.X(I2).EQ.ANUM.OR.X(I3).EQ.ANUM) THEN
         CALL CPRINT(
     +    'Error in CARTX2 -- One of the antecedent atoms is undefined.'
     +    )
 
         WRITE (BUFFER,9001) I1, X(I1), Y(I1), Z(I1)
         CALL CPRINT(BUFFER)
 
         WRITE (BUFFER,9002) I2, X(I2), Y(I2), Z(I2)
         CALL CPRINT(BUFFER)
 
         WRITE (BUFFER,9003) I3, X(I3), Y(I3), Z(I3)
         CALL CPRINT(BUFFER)
 
         CALL DIE
      ENDIF
C
C     Generate position of X4 with everything else easily lined up.
C
      STHT = SIN(PI-THETA)
      CTHT = COS(PI-THETA)
      SPHI = SIN(PHI)
      CPHI = COS(PHI)
      BSIN = BOND*STHT
      X4 = BOND*CTHT
      Y4 = BSIN*CPHI
      Z4 = BSIN*SPHI
C
C     Translate X1 and X2 so that X3 would be at the origin.
C
      X3 = X(I3)
      Y3 = Y(I3)
      Z3 = Z(I3)
      X1 = X(I1) - X3
      Y1 = Y(I1) - Y3
      Z1 = Z(I1) - Z3
      X2 = X(I2) - X3
      Y2 = Y(I2) - Y3
      Z2 = Z(I2) - Z3
C
C     Rotate X1 by rotation of X2 to the origin.
C
      LXZ22 = X2*X2 + Z2*Z2
      L2 = SQRT(LXZ22+Y2*Y2)
      LXZ2 = SQRT(LXZ22)
      IF (L2.LT.ETA) THEN
         CALL CPRINT(
     +'0***** Warning from CARTX ***** Atom 2 and Atom 3 are too close.'
     +)
         OVL2 = 1.0/ETA
      ELSE
         OVL2 = 1.0/L2
      ENDIF
      IF (LXZ2.LT.ETA) THEN
         XX1 = X1
         X2O = 1.0
         Z2O = 0.0
      ELSE
         OVLXZ2 = 1.0/LXZ2
         X2O = X2*OVLXZ2
         Z2O = Z2*OVLXZ2
         XX1 = X1*X2O + Z1*Z2O
         Z1 = Z1*X2O - X1*Z2O
      ENDIF
      XZ2O = LXZ2*OVL2
      Y2O = Y2*OVL2
      X1 = -XX1*XZ2O - Y1*Y2O
      Y1 = XX1*Y2O - Y1*XZ2O
C
C     Rotate X4 by inverse of rotation which takes the transformed X1 to
C     the XY plane by rotation about the x axis.
C
      LYZ1 = SQRT(Y1*Y1+Z1*Z1)
      OVLYZ1 = 1.0/LYZ1
      Y1O = Y1*OVLYZ1
      Z1O = Z1*OVLYZ1
      YY4 = Y1O*Y4 - Z1O*Z4
      ZZ4 = Y1O*Z4 + Z1O*Y4
C
C     Rotate X4 by inverse of X2 rotation to -X axis.
C
      XX4 = Y2O*YY4 - XZ2O*X4
      Y4 = -XZ2O*YY4 - Y2O*X4
      X4 = X2O*XX4 - Z2O*ZZ4
      Z4 = Z2O*XX4 + X2O*ZZ4
      X(I4) = X4 + X3
      Y(I4) = Y4 + Y3
      Z(I4) = Z4 + Z3
C
      RETURN
C
 9001 FORMAT (' Atom 1 index = ',I5,3F10.3)
 9002 FORMAT (' Atom 2 index = ',I5,3F10.3)
 9003 FORMAT (' Atom 3 index = ',I5,3F10.3)
C
      END

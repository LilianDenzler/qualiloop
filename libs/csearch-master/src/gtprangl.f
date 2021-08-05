CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION GTPRANGL(IT,JT,KT,PARM_NO,IAC,KCT,NCT,NATC,TYPE,
     +                          CTB)
C
C     Returns the bond angle for the angle between atoms IT, JT, & KT.
C
      IMPLICIT NONE
 
 
      INTEGER IT, JT, KT, I, ICTB
      INTEGER*2 PARM_NO(100,100), IAC(*)
      INTEGER KCT(*), NCT, NATC, TYPE(*)
      REAL CTB(*)
      INTEGER J, K, CODE, NINDX
 
      CHARACTER*100 BUFFER
C
      IF (IT.EQ.0.OR.JT.EQ.0.OR.KT.EQ.0) THEN
         WRITE (BUFFER,9001)
         CALL CPRINT(BUFFER)
         CALL DIE
      ENDIF
      I = IAC(IT)
      J = IAC(JT)
      K = IAC(KT)
      CODE = PARM_NO(I,K)
      CODE = NATC*CODE + J
      ICTB = NINDX(CODE,KCT,NCT)
      IF (ICTB.EQ.0) THEN
         WRITE (BUFFER,9002) TYPE(IT), TYPE(JT), TYPE(KT)
         CALL CPRINT(BUFFER)
         CALL DIE
      ENDIF
      GTPRANGL = CTB(ICTB)
      RETURN
 9001 FORMAT ('Error in GTPRANGL: Atom index is zero')
 9002 FORMAT ('0Error in GTPRANGL -- Unable to find parameters for',
     +        ' angle between ',A4,' and ',A4,' and ',A4)
      END

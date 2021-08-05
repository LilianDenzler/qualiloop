C     Copyright (c) 1987 Robert E. Bruccoleri
C     Copying of this software, in whole or in part, is permitted
C     provided that the copies are not made for commercial purposes,
C     appropriate credit for the use of the software is given, this
C     copyright notice appears, and notice is given that the copying
C     is by permission of Robert E. Bruccoleri. Any other copying
C     requires specific permission.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION GETPARBOND(IB,JB,PARM_NO,IAC,KCB,NCB,TYPE,CBB)
C
C     Returns the bond length for the bond between atoms IB and JB.
C
      IMPLICIT NONE
 
      INTEGER IB, JB, I, ICBB
      INTEGER*2 PARM_NO(100,100), IAC(*)
      INTEGER KCB(*), NCB, TYPE(*)
      REAL CBB(*)
      INTEGER J, CODE, NINDX
 
      CHARACTER*100 BUFFER
C
      IF (IB.EQ.0.OR.JB.EQ.0) THEN
         CALL CPRINT('Error in GETPARBOND: Atom index is zero')
         CALL DIE
      ENDIF
      I = IAC(IB)
      J = IAC(JB)
      CODE = PARM_NO(I,J)
      ICBB = NINDX(CODE,KCB,NCB)
      IF (ICBB.EQ.0) THEN
         WRITE (BUFFER,9001) TYPE(IB), TYPE(JB)
         CALL CPRINT(BUFFER)
         CALL DIE
      ENDIF
      GETPARBOND = CBB(ICBB)
      RETURN
 9001 FORMAT ('Error in GETPARBOND -- Unable to find parameters for',
     +        ' bond between ',A4,' and ',A4)
      END

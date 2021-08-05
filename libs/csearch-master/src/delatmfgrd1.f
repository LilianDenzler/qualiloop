CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE delatmfgrd1(IND,SPACE_GRID,INGRID,NEXTHD,
     +                                  CLSHD,CLSTL,NEXTCLS,CLSATM)
C
C     Deletes atom IND from the grid.
C
      include "impnone.inc"
      include "grid.inc"
      INTEGER IND
      INTEGER*2 SPACE_GRID(NGRIDX,NGRIDY,NGRIDZ)
      LOGICAL INGRID(*)
      INTEGER NEXTHD(*), CLSHD(*), CLSTL(*), NEXTCLS(*), CLSATM(*)
      CHARACTER*512 BUFFER
      include "params.inc"
      include "coords.inc"
      include "dbg.inc"
C
      INTEGER IX, IY, IZ, ISPACE
      LOGICAL OUTOFBOUND
C
      IF (.NOT.INGRID(IND)) THEN
         PRINT *, 'Using format 9001'
         WRITE (BUFFER,9001) IND
         CALL CPRINT(BUFFER)
         CALL DIE
      ENDIF
      INGRID(IND) = .FALSE.
      CALL COMPSPCGRD(IND,IX,IY,IZ,OUTOFBOUND)
      IF (OUTOFBOUND) THEN
         PRINT *, 'Using format 9002'
         WRITE (BUFFER,9002) IX, IY, IZ
         CALL CPRINT(BUFFER)
         CALL DIE
      ENDIF
      ISPACE = SPACE_GRID(IX,IY,IZ)
      IF (ISPACE.EQ.0) THEN
         PRINT *, 'Using format 9003'
         WRITE (BUFFER,9003) IND
         CALL CPRINT(BUFFER)
         CALL DIE
      ENDIF
      CALL SCHDLS(CLSHD(ISPACE),CLSTL(ISPACE),FREECLS,NEXTCLS,CLSATM,
     +            IND)
      IF (CLSHD(ISPACE).EQ.0) THEN
         SPACE_GRID(IX,IY,IZ) = 0
         NEXTHD(ISPACE) = FREEHD
         FREEHD = ISPACE
      ENDIF
      RETURN
 9001 FORMAT ('Error in delatmfgrd1 -- Atom',I6,
     +        ' deleted twice from CLOSE CONTACTS.')
 9002 FORMAT ('Error in delatmfgrd1 --',
     +        ' Bad indices for SPACE_GRID = ',3I12)
 9003 FORMAT ('Error in delatmfgrd1 -- Atom',I12,
     +        ' not found in SPACE_GRID.')
      END

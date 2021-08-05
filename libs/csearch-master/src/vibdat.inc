**!!VIBDAT.INS
C
C
C -- Modified:-
C      13-Sep-1990      JPH    'temp' renamed to 'vtemp' - stops
C                              conflict with other 'temp' in common
C                              block 'mdopt'.
C
C -- ????
      DOUBLE PRECISION VTEMP
      COMMON /VIBDAT/  VTEMP
C
C -- ????
      REAL SYMG
C -- ????
      REAL SYM
      COMMON /VIBDAT/ SYMG , SYM ( 160 )
C
      SAVE /VIBDAT/
C

**!!ZMOLIN.INS
C
C -- ????
      DOUBLE PRECISION HBDISA
C -- ????
      DOUBLE PRECISION UCHRG
      COMMON /ZMOLIN/  HBDISA , UCHRG
C
C -- ????
      INTEGER ICART
C -- ????
      INTEGER ICONCT
C -- ????
      INTEGER IPRLIB
C -- ????
      INTEGER IRSLIB
C -- ????
      INTEGER IUPCON
C -- ????
      INTEGER IUPCOR
C -- ????
      INTEGER KFIX
C -- ????
      INTEGER KGRPMN
C -- ????
      INTEGER KINPPR
C -- ????
      INTEGER KOPLNG
C -- ????
      INTEGER KTRNLB
C -- ????
      INTEGER LLSQ
C -- ????
      INTEGER NVSET
C -- ????
      INTEGER ITLRMS ( 20 )
C -- ????
      INTEGER IVSET ( MXVSET )
C -- ????
      INTEGER LINT ( MXLINT )
C -- ????
      INTEGER LVSET ( MXVSET )
C -- Cell data read ?
      INTEGER ICELLD
C
      COMMON /ZMOLIN/  ICART  , ICONCT , IPRLIB , IRSLIB , IUPCON ,
     +                 IUPCOR , KFIX   , KGRPMN , KINPPR , KOPLNG ,
     +                 KTRNLB , LLSQ   , NVSET  , ITLRMS , IVSET  ,
     +                 LINT   , LVSET  , ICELLD
C
      SAVE /ZMOLIN/
C

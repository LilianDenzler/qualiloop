**!!MDTEMP.INS
C
C -- Store the atomic temperatures during dynamics.
C
      DOUBLE PRECISION ATMTMP ( MXATM )
      COMMON /MDTEMP/ ATMTMP
C
      SAVE /MDTEMP/

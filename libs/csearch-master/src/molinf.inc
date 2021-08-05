C
C -- Name:-
C	MOLINF.INC
C
C -- Function:-
C	Molcule information for eureka.
C
C
C -- Modified:-
C	J.P.Holland	16-JUL-1991	Oxford Molecular Ltd.
C
C
C
      INTEGER*2 MOLXYZ, NATMOL, MOLRES, NRSMOL, MOLGRP, NGPMOL, MOLBP,
     +          NBMOL,  MOLTP,  NTMOL,  MOLPP,  NPMOL,  MOLOP,  NOPMOL,
     +          MOLTTP, NTTMOL, MOLBND, NBDMOL
C
      COMMON /MOLPTR/ MOLXYZ(MXMOL), NATMOL(MXMOL), MOLRES(MXMOL),
     +                NRSMOL(MXMOL), MOLGRP(MXMOL), NGPMOL(MXMOL),
     +                MOLBP(MXMOL),  NBMOL(MXMOL),  MOLTP(MXMOL),
     +                NTMOL(MXMOL),  MOLPP(MXMOL),  NPMOL(MXMOL),
     +                MOLOP(MXMOL),  NOPMOL(MXMOL), MOLTTP(MXMOL),
     +                NTTMOL(MXMOL), MOLBND(MXMOL), NBDMOL(MXMOL)
C
      SAVE /MOLPTR/

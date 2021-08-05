**!!MNOPT.INS
C
C -- Number of cycles of STEEPEST DESCENT minimiser
      INTEGER NSTEEP
C -- Number of cycles of CONJUGATE GRADIENT minimiser
      INTEGER NCJGR
C -- Number of cycles of QMIN ( VA09A ) minimiser
      INTEGER NMINO
C -- Number of cycles of NEWTON-RAPHSONminimiser
      INTEGER NRAF
C -- Number of iterations between saving structure
      INTEGER NSAVE
C
C
C -- VA09A mode
      INTEGER IVMODE
C -- VA09A print options - print infornmation every 'IVPR' iterations
      INTEGER IVPR
C -- VA09A restart option
C      = 1  write out estimated 2nd deriv array at end of minimisation
C        2  read in at start and write out at end
      INTEGER KVARST
C
C -- NEWTON-RAPHSON conjugate gradient options
C      = 1   self-correcting conj gradient minimisation
C      else  HSCG
      INTEGER LSLFCR
C -- NEWTON-RAPHSON max number of interpolations in line search
C      - *** NOT USED ***
      INTEGER INTMXA
C -- NEWTON-RAPHSON restart every 'IRESTR' cycles
      INTEGER IRESTR
C
C
C -- CONJUGATE GRADIENT print option - print every 'IFRQPR' iterations
      INTEGER IFRQPR
C -- CONJUGATE GRADIENT update formula used to update conj. grad. factor
C      = 1  Fletcher-Reeves
C        2  Polak-Ribiere
C        3  Hestenes-Stiefel
      INTEGER ISWTCH
C -- Total number of iteration in minimisers
      INTEGER NCOUNT
      COMMON /MNOPT/ NSTEEP , NCJGR  , NMINO  , NRAF   , NSAVE  ,
     +               IVMODE , IVPR   , LSLFCR , INTMXA , IRESTR ,
     +               IFRQPR , ISWTCH , KVARST , NCOUNT
C
C -- STEEPEST DESCENT scaling factor for minimization step
      DOUBLE PRECISION STEP0
C -- STEEPEST DESCENT scaling factor for step size if energy decreased
C      - step size is enlarged if energy decreased last step
      DOUBLE PRECISION FCTINC
C -- STEEPEST DESCENT scaling factor for step size if energy increased
C      - step size is reduced if energy increased last step
      DOUBLE PRECISION FCTDEC
C -- STEEPEST DESCENT derivative termination option
C      - stop minimization when all derivatives less than 'DERSTP'
      DOUBLE PRECISION DERSTP
C -- STEEPEST DESCENT rms termination option
C      - stop minimization when RMS is less than 'RMSSTP'
      DOUBLE PRECISION RMSSTP
C
C
C -- QMIN initial scaling factor
      DOUBLE PRECISION DFN
C -- QMIN termination 1 - finish when all derivs less than 'DERMAX'
      DOUBLE PRECISION DERMAX
C -- QMIN termination 2 - finish when RMS less than 'RMSLMT'
      DOUBLE PRECISION RMSLMT
C -- QMIN initial value of 2nd derivative matrix diagonal
      DOUBLE PRECISION HSCALE
C -- QMIN termination 3 - finish when change in step less than 'EXSMST'
      DOUBLE PRECISION EXSMST
C
C
C -- NEWTON-RAPHSON minimum fraction of old fn. value decrease expected
C      - *** NOT USED ***
      DOUBLE PRECISION GESFCA
C -- NEWTON-RAPHSON min allowed line search step size
C      - *** NOT USED ***
      DOUBLE PRECISION STPMNA
C -- NEWTON-RAPHSON max allowed line search step size
C      - *** NOT USED ***
      DOUBLE PRECISION STPMXA
C -- NEWTON-RAPHSON accuracy of line search control ( 0.0 to 0.1 )
C      - *** NOT USED ***
      DOUBLE PRECISION CAPPAA
C -- NEWTON-RAPHSON termination if rel. fn. increase exceeds 'DELTAA'
C      - *** NOT USED ***
      DOUBLE PRECISION DELTAA
C -- NEWTON-RAPHSON termination if tot coord change less than 'EPSLNA'
C      - *** NOT USED ***
      DOUBLE PRECISION EPSLNA
C -- NEWTON-RAPHSON steepest descent while fn. decrease exceeds 'FASTA'
C      - *** NOT USED ***
      DOUBLE PRECISION FASTA
C -- NEWTON-RAPHSON restart if rel. fn. decrease less than 'SLOWA'
C      - *** NOT USED ***
      DOUBLE PRECISION SLOWA
C -- NEWTON-RAPHSON restart if angle between grad. and conj. grad
C      - *** NOT USED ***
C      exceeds 'ANGMXA'
      DOUBLE PRECISION ANGMXA
C -- NEWTON-RAPHSON exit if fn. value less than 'FCTMNA'
C      - *** NOT USED ***
      DOUBLE PRECISION FCTMNA
C -- NEWTON-RAPHSON exit if rel energy change less than 'ACCNR'
      DOUBLE PRECISION ACCNR
C -- NEWTON-RAPHSON ignore eigenvalues less than 'EXEIG'
      DOUBLE PRECISION EXEIG
C -- NEWTON-RAPHSON scaling factor for minimization step
      DOUBLE PRECISION STPNEW
C
C
C -- CONJUGATE GRADIENT accuracy of line search control ( 0.0 TO 0.1 )
      DOUBLE PRECISION SPREC
C -- CONJUGATE GRADIENT upper bound on multiplicative increase in
C      search factor
      DOUBLE PRECISION EXTBND
C -- CONJUGATE GRADIENT termination when sum of squares of derivs less
C      than 'TOL'
      DOUBLE PRECISION TOL
C -- CONJUGATE GRADIENT estimate of expected reduction in energy
      DOUBLE PRECISION RFN
C -- CONJUGATE GRADIENT terminate when all derivs less than 'DERCNJ'
      DOUBLE PRECISION DERCNJ
C -- CONJUGATE GRADIENT terminate when RMS less than 'RMSCNJ'
      DOUBLE PRECISION RMSCNJ
C
C
      COMMON /MNOPT/   STEP0  , FCTINC , FCTDEC , DFN    , DERMAX ,
     +                 HSCALE , EXSMST , GESFCA , STPMNA , STPMXA ,
     +                 CAPPAA , DELTAA , EPSLNA , FASTA  , SLOWA  ,
     +                 ANGMXA , FCTMNA , ACCNR  , EXEIG  , SPREC  ,
     +                 EXTBND , TOL    , RFN    , RMSLMT , DERSTP ,
     +                 RMSSTP , RMSCNJ , DERCNJ , STPNEW
C
      SAVE /MNOPT/
C

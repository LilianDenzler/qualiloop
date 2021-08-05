**!!MDOPT.INS
C
C -- Read unit number.
      INTEGER IREAD
C -- Print unit number.
      INTEGER IPRINT
C -- Screen unit number
      INTEGER ISCRN
C -- Restart file unit number.
      INTEGER LINRST
C -- ????
      INTEGER LINITV
C -- ????
      INTEGER LOPRST
C -- History file unit number.
      INTEGER LHIS
C -- Print control.
C      = 0  Print initial and final energies.
C        1  Plus coordinates, velocities and accelerations per molecule.
C        2  Plus coordinates, velocities and accelerations per atom.
      INTEGER KPRINT
C -- Averages print control.
C      = 0  No averages printout.
C        1  Total averages printed every 'NSTPAV' steps.
C        2  Plus averages of all energy components every 'NSTPAV' steps.
      INTEGER KAVPR
C -- Step print controls.
C      = 0  No printout every step.
C        1  Total energies printed every step.
C        2  Plus energy components.
C        3  Plus linear/angular velocity.
C        4  Plus coordinates, velocities, accelerations per molecule.
C        5  Plus coordinates, velocities, accelerations per atom.
C        6  Plus coords. etc printed in Gear algorithm for debugging.
      INTEGER KSTPPR
C -- Restart run ?
C      = 0  Restart run
C        1  Initial run
      INTEGER KINIT
C -- Compute average energies and temperatures every 'NSTPAV' steps.
      INTEGER NSTPAV
C -- Number of steps to perform.
      INTEGER NSTEPS
C -- Write history every 'NWRSTP' steps.
      INTEGER NWRSTP
C -- Clear higher order derivatives on rescale ( Gear only ).
      INTEGER KCLRH
C -- Velocity rescale options.
C      = 0  No rescaling.
C        1  Temperature rescaling.
C        2  Total energy rescaling.
C        3  Single rescale at start of run.
      INTEGER KVRSCL
C -- Rescale "hot" atoms ?
C      = 0  No
C        1  Rescale atoms with temp greater than 'ATTEMP'.
      INTEGER KINDV
C -- Exact rescale ?
C      = 0  Rescale velocities by 1.0+/-VSCALE only.
C        1  Exact rescale.
      INTEGER KEXACT
C -- Test temperature each step for rescale ?
C      = 0  Yes
C        1  Test every 'NSTPAV' steps.
      INTEGER KRSCLV
C -- Restart option - time step rescale option.
C      = 1  Change time step to new 'TIMSTP' value.
      INTEGER KTSSCL
C -- Integration algorithm.
C      = 1  Gear.
C        2  Verlet.
C        3  Beeman.
C        4  Leapfrog.
      INTEGER KALG
C -- ????
      INTEGER NORDRI
C -- Order for Gear algorithm.
      INTEGER NORDER
C -- Number of corrector steps for Gear algorithm.
      INTEGER NCORCT
C -- Initial setup file unit number ( not used ).
      INTEGER KINITV
C -- ????
      INTEGER KRSTFC
C -- Equilibrate temp at each atom individually by velocity rescaling ?
C      = 0  No.
C        1  Yes.
      INTEGER KEQUIL
C -- Constant temperature bath ?
C      = 0  No.
C        1  Yes.
      INTEGER KTEMPB
C -- Constant pressure bath ?
C      = 0  No.
C        1  Yes.
      INTEGER KPRESB
C -- Allow shape of unit cell to change when scaling coords for const P.
C      = 0  No  - only unit cell lengths are changed.
C        1  Yes - unit cell lengths and angles may be changed.
      INTEGER KSHAPE
C -- Internal pressure calculated using virial based on :
C      = 0  centre of mass of each molecule.
C        1  coordinates of each atom.
      INTEGER KATPRS
C -- Pressure scaling factor calculated using virial based on :
C      = 0  centre of mass of each molecule.
C        1  coordinates of each atom.
      INTEGER KATPFC
C -- Constant pressure scaling applied to :
C      = 0  centre of mass of each molecule as a rigid body.
C        1  coordinates of each atom.
      INTEGER KATSCL
C -- Bond lengths constrained to equilibrium lengths ?
C      = 0  No.
C        1  Yes.
      INTEGER KCNSTR
C -- Write molecular dynamics history file ?
C      = 0  Yes.
C        1  No.
      INTEGER KMDHIS
      COMMON /MDOPT/ IREAD  , IPRINT , LINRST , LINITV , LOPRST , LHIS ,
     +               KPRINT , KAVPR  , KSTPPR , KINIT  , NSTPAV ,
     +               NSTEPS , NWRSTP , KCLRH  , KVRSCL , KINDV  ,
     +               KEXACT , KRSCLV , KTSSCL , KALG   , NORDRI ,
     +               NORDER , NCORCT , KINITV , KRSTFC , KEQUIL ,
     +               KTEMPB , KPRESB , KSHAPE , KATPRS , KATPFC ,
     +               KATSCL , KCNSTR , KMDHIS , ISCRN
C
C -- ????
*-IF FPS
*-      WORD ITITLR , JTITLR
*-ELSE
      INTEGER ITITLR , JTITLR
*-ENDIF
      COMMON /MDOPT/ ITITLR ( 20 ) , JTITLR ( 20 )
C
C -- Stop if within 'FINISH' seconds of max CPU time ( not used ).
      DOUBLE PRECISION FINISH
C -- Stop if total energy changes by more than 'DEMAX'.
      DOUBLE PRECISION DEMAX
C -- Energy rescale - total energy required.
      DOUBLE PRECISION TOTEN
C -- Energy rescale - rescale when energy changes by more than 'DTOTEN'.
      DOUBLE PRECISION DTOTEN
C -- Temp rescale - constant temperature required.
      DOUBLE PRECISION TEMP
C -- Temp rescale - rescale when temp changes by more than 'DTEMP'.
      DOUBLE PRECISION DTEMP
C -- Temp rescale - scale factor for 'KEXACT' = 0.
C    Single rescale - velocities multiplied directly by this factor.
      DOUBLE PRECISION VSCALE
C -- Temp rescale - cutoff temp for "hot" atoms.
      DOUBLE PRECISION ATTEMP
C -- Temp rescale - scale factor to reduce velocities of "hot" atoms.
      DOUBLE PRECISION ATSCAL
C -- Input time step
      DOUBLE PRECISION TMSTPI
C -- Old time step ( if restart run ).
      DOUBLE PRECISION TMSTPO
C -- Version number.
      DOUBLE PRECISION VERSHN
C -- Seed for random velocities.
      DOUBLE PRECISION DSEED
C -- Temperature for initial random velocities
      DOUBLE PRECISION TEMP0
C -- Temperature required for const T bath.
      DOUBLE PRECISION TEMPB
C -- Coupling timestep for const T bath.
      DOUBLE PRECISION TIMTMP
C -- Pressure for const P bath.
      DOUBLE PRECISION PRESSB
C -- Coupling timestep for const P bath.
      DOUBLE PRECISION TIMPRS
C -- Isothermal compressibility for const P bath.
      DOUBLE PRECISION CMPRSB
      COMMON /MDOPT/   FINISH , DEMAX  , TOTEN  , DTOTEN , TEMP   ,
     +                 DTEMP  , VSCALE , ATTEMP , ATSCAL , TMSTPI ,
     +                 TMSTPO , VERSHN , DSEED  , TEMP0  , TEMPB  ,
     +                 TIMTMP , PRESSB , TIMPRS , CMPRSB
C
C
C -- Pre-warminig run before main run ?
      INTEGER KINITR
C -- Temperature to pre-warm to.
      DOUBLE PRECISION TEMPR
C -- Max steps allowed for pre-warming run.
      INTEGER NSTEPR
      INTEGER ISPARE1
C -- Write restart file every n steps.
      INTEGER NWRSTT
C -- Average temperature
      DOUBLE PRECISION AVTEMP
C -- Total energy
      DOUBLE PRECISION AVTOT
C -- Kinetic energy
      DOUBLE PRECISION AVTK
C -- Potential energy
      DOUBLE PRECISION AVTP
      COMMON /MDOPT/ KINITR , TEMPR , NSTEPR , NWRSTT ,
     +               AVTEMP , AVTOT , AVTK  , AVTP
C
      SAVE /MDOPT/

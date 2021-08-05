**!!FFOPT.INS
C
C
C -- Use Morse or harmonic potential for bond strain
C      = 0  harmonic
C        1  Morse     (default)
      INTEGER IMORSE
C -- ** NOT USED ** Use Urey-Bradley force field
C      = 0  No        (default)
C        1  Yes
      INTEGER IUB
C -- ** NOT USED ** Use valence force field
C      = 0  No
C        1  Yes       (default)
      INTEGER IVFF
C -- ** NOT USED ** Use Allinger force field
C      = 0  No        (default)
C        1  Yes
      INTEGER IALNGR
C -- Exclude 1-4 interactions
C      = 0  No
C        1  Yes       (default)
      INTEGER IEXC
C -- Exclude cross-terms
C      = 0  No        (default)
C        1  Yes
      INTEGER ICROSS
C -- Use 10-12 potential for hydrogen bonds ( Allinger only )
C      = 0  No        (default)
C        1  Yes
      INTEGER K1012
C -- Use energy cutoff ?
      INTEGER ICUT
      COMMON /FFOPT/ IMORSE , IUB , IVFF , IALNGR , IEXC , ICROSS ,
     +               K1012  , ICUT
C
C -- Include all residues within 'CUTOFF' in neighbour list
      DOUBLE PRECISION CUTOFF
C -- ** NOT USED ** Include all interactions within 'CUTREP' in repulsn energy
      DOUBLE PRECISION CUTREP
C -- ** NOT USED ** Include all interactions within 'CUTDSP' in dispersn energy
      DOUBLE PRECISION CUTDSP
C -- ** NOT USED ** Include all interactions within 'CUTCOU' in charge energy
      DOUBLE PRECISION CUTCOU
C -- Include all interactions within 'CUTDIS' in energy
      DOUBLE PRECISION CUTDIS
C -- Switching fn. reduces non-bond energy to zero outside 'SWTDIS'
      DOUBLE PRECISION SWTDIS
C -- Scaling factor for 1-4 electrostatic interactions
      DOUBLE PRECISION SF14QQ
C -- Scaling factor for 1-4 non-bond interactions
      DOUBLE PRECISION SF14LJ
      COMMON /FFOPT/   CUTOFF , CUTREP , CUTDSP , CUTCOU , CUTDIS ,
     +                 SWTDIS , SF14QQ , SF14LJ
C
      SAVE /FFOPT/
C

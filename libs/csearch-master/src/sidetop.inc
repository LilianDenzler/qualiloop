C***********************************************************************
C*      NAME: SIDETOP                                                  *
C*  FUNCTION: To declare variables relating to sidechain topology      *
C* COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
C*---------------------------------------------------------------------*
C*    AUTHOR: Robert Williams                                          *
C*      DATE: 05/10/92                                                 *
C*---------------------------------------------------------------------*
C*    INPUTS:                                                          *
C*   OUTPUTS:                                                          *
C*    LOCALS:                                                          *
C*   GLOBALS:                                                          *
C*     CALLS:                                                          *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DD/MM/YY   INITS   COMMENTS                                         *
C***********************************************************************
C
C   Variable     Array bounds            Description
C   --------     ------------            -----------

C   IADD              -       Flag for additions to free torsion
C   IFIXED            -       Flag indicating fixed torsion
C   IFREE             -       Torsion flag indicating independent atom
C   NSCATM            -       The number of sidechain atoms
C   NSCLMP            -       The number of sidechain clumps
C   NSCRES            -       The number of sidechain residues

C   SC_ATOM_BLD    MXSCAT     IUPAC name of atom to be constructed
C   SC_ANTE1_BLD   MXSCAT     Neighbour IUPAC named atom for building
C   SC_ANTE2_BLD   MXSCAT     Next adjacent IUPAC atom for building
C   SC_ANTE3_BLD   MXSCAT     Third adjacent IUPAC atom for building
C   SC_BOND_BLD    MXSCAT     Bond for building an atom
C   SC_ANGLE_BLD   MXSCAT     Construction angle for each atom
C   SC_CODE_BLD    MXSCAT     Flag describing the use of the torsion angle
C   SC_TORS_BLD    MXSCAT     Torsion angle for building (FREE=> dont use)
C   SC_ATOM_PART   MXCLMP     Partition of atoms into clumps
C   SC_CLUMP_PART  MXSRES     Partition of clumps into residues 
C   SC_FREE_ATOM   MXCLMP     The free atom in the clump
C   SC_OFFSET      MXSCAT     Offset angle for torsion ADD procedure
C   SC_RESNAME     MXSRES     The residue names
C   SC_SPECIAL     MXSRES     TRUE if consistency checks are ignored
C   SC_SYMMETRY    MXCLMP     The Rotational symmetry of each clump

C Where:
C
C   MXSRES                 Maximum number of sidechain residues
C   MXCLMP                 Maximum number of sidechain clumps
C   MXSCAT                 Maximum number of sidechain atoms

      INTEGER MXSRES, MXCLMP, MXSCAT
      PARAMETER (MXSRES = 30, MXCLMP = 75, MXSCAT = 200)

      INTEGER IADD, IFIXED, IFREE, NSCATM, NSCLMP, NSCRES
      DATA IADD, IFIXED, IFREE /1,2,3/

      INTEGER SC_ATOM_BLD(MXSCAT), SC_ANTE1_BLD(MXSCAT),
     +        SC_ANTE2_BLD(MXSCAT), SC_ANTE3_BLD(MXSCAT),
     +        SC_CODE_BLD(MXSCAT),  SC_ATOM_PART(MXCLMP+1), 
     +        SC_CLUMP_PART(MXSRES+1), SC_FREE_ATOM(MXCLMP),
     +        SC_RESNAME(MXSRES), SC_SPECIAL(MXSRES), 
     +        SC_SYMMETRY(MXCLMP) 

      REAL SC_ANGLE_BLD(MXSCAT), SC_BOND_BLD(MXSCAT),
     +     SC_TORS_BLD(MXSCAT), SC_OFFSET(MXSCAT)

      COMMON /SIDETOP/ SC_BOND_BLD, SC_ANGLE_BLD, SC_TORS_BLD,
     +                 SC_OFFSET, NSCATM, NSCLMP, NSCRES, SC_CODE_BLD,
     +                 SC_ATOM_PART, SC_CLUMP_PART, SC_SPECIAL,
     +                 SC_SYMMETRY, SC_ANTE1_BLD, SC_ANTE2_BLD,
     +                 SC_ANTE3_BLD, SC_RESNAME, SC_FREE_ATOM,
     +                 SC_ATOM_BLD
      SAVE  /SIDETOP/

/***********************************************************************
 *      NAME: VALUES                                                   *
 *  FUNCTION: To declare the most frequently used system values        *
 * COPYRIGHT: (C) OXFORD MOLECULAR LIMITED, 1992                       *
 *---------------------------------------------------------------------*
 *    AUTHOR: Robert Williams                                          *
 *      DATE: 05/10/92                                                 *
 *---------------------------------------------------------------------*
 *    INPUTS:                                                          *
 *   OUTPUTS:                                                          *
 *    LOCALS:                                                          *
 *   GLOBALS:                                                          *
 *     CALLS:                                                          *
 *---------------------------------------------------------------------*
 * MODIFICATION RECORD                                                 *
 * DD/MM/YY   INITS   COMMENTS                                         *
 ***********************************************************************
     Variable              Description
     --------              -----------
 
      NATOMS        The number of atoms in the system
      NRES          The number of residues
      NSEGS         The number of segments
      NBONDS        The number of bonds in the system
      NANGS         The number of bond angles in the system
      NPTORS        The number of proper torsion angles
      NITORS        The number of improper torsion angles
      NHBS          The number of hydrogen bonds
      NNBS          The number of non-bond pairs
      NDONAT        The number of hydrogen bond donor atoms
      NACCAT        The number of hydrogen bond acceptor atoms
      NBPAR         The number of bond parameters
      NAPAR         The number of angle parameters
      NPTPAR        The number of proper torsion parameters
      NITPAR        The number of improper torsion parameters
      NHBPAR        The number of hydrogen parameters
      NATYPS        The number of atom types
      NBAUTO        Flag to generate non-bonded exclusions
*/
 
#ifdef __MAIN__
struct {

      int natoms, nres, nsegs, nbonds, nangs, nptors, nitors, nhbs,
          nnbs, ndonat, naccat, nbpar, napar, nptpar, nitpar,
          nhbpar, natyps, nbauto;

      } values;
#else
extern struct {

      int natoms, nres, nsegs, nbonds, nangs, nptors, nitors, nhbs,
          nnbs, ndonat, naccat, nbpar, napar, nptpar, nitpar,
          nhbpar, natyps, nbauto;

      } values;

#endif


#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif
 
/* This routine removes all internal coordinates which have atoms outside
   the range of the atoms used in the current segment. The current segment
   is defined by nictot and nseg. Deletion is done by assigning all the
   internal coordinates an atom index of 0.

   04.08.92 Some rewriting.   By:   ACRM.
*/
 
void trimic(
void
)
{
   int   FirstAtom,
         FirstBond,
         FirstAngle,
         FirstProper,
         FirstImproper,
         FirstNBond,
         InSegment,
         i;
 
   FirstAtom     = pstruct.segndx[values.nsegs-1][1] + 1;
   FirstBond     = pstruct.segndx[values.nsegs-1][2] + 1;
   FirstAngle    = pstruct.segndx[values.nsegs-1][3] + 1;
   FirstProper   = pstruct.segndx[values.nsegs-1][4] + 1;
   FirstImproper = pstruct.segndx[values.nsegs-1][5] + 1;
   FirstNBond    = pstruct.segndx[values.nsegs-1][6] + 1;
 
   for(i=FirstBond; i<=values.nbonds; i++)
   {
      InSegment = (pstruct.atbnd1[i-1] >= FirstAtom &&
                   pstruct.atbnd1[i-1] <= values.natoms &&
                   pstruct.atbnd2[i-1] >= FirstAtom &&
                   pstruct.atbnd2[i-1] <= values.natoms) ? 
                   TRUE : FALSE;
      if(!InSegment)
      {
         pstruct.atbnd1[i-1] = 0;
         pstruct.atbnd2[i-1] = 0;
      }
   }
 
   for(i=FirstAngle; i<=values.nangs; i++)
   {
      InSegment = (pstruct.atang1[i-1] >= FirstAtom &&
                   pstruct.atang1[i-1] <= values.natoms &&
                   pstruct.atang2[i-1] >= FirstAtom &&
                   pstruct.atang2[i-1] <= values.natoms &&
                   pstruct.atang3[i-1] >= FirstAtom &&
                   pstruct.atang3[i-1] <= values.natoms) ?
                   TRUE : FALSE;
      if(!InSegment)
      {
         pstruct.atang1[i-1] = 0;
         pstruct.atang2[i-1] = 0;
         pstruct.atang3[i-1] = 0;
      }
   }
 
   for(i=FirstProper; i<=values.nptors; i++)
   {
      InSegment = (pstruct.attor1[i-1] >= FirstAtom &&
                   pstruct.attor1[i-1] <= values.natoms &&
                   pstruct.attor2[i-1] >= FirstAtom &&
                   pstruct.attor2[i-1] <= values.natoms &&
                   pstruct.attor3[i-1] >= FirstAtom &&
                   pstruct.attor3[i-1] <= values.natoms &&
                   pstruct.attor4[i-1] >= FirstAtom &&
                   pstruct.attor4[i-1] <= values.natoms ) ? 
                   TRUE : FALSE;
      if(!InSegment)
      {
         pstruct.attor1[i-1] = 0;
         pstruct.attor2[i-1] = 0;
         pstruct.attor3[i-1] = 0;
         pstruct.attor4[i-1] = 0;
      }
   }
 
   for(i=FirstImproper; i<=values.nitors; i++)
   {
      InSegment = (pstruct.atimp1[i-1] >= FirstAtom &&
                   pstruct.atimp1[i-1] <= values.natoms &&
                   pstruct.atimp2[i-1] >= FirstAtom &&
                   pstruct.atimp2[i-1] <= values.natoms &&
                   pstruct.atimp3[i-1] >= FirstAtom && 
                   pstruct.atimp3[i-1] <= values.natoms &&
                   pstruct.atimp4[i-1] >= FirstAtom &&
                   pstruct.atimp4[i-1] <= values.natoms ) ? 
                   TRUE : FALSE;
      if(!InSegment)
      {
         pstruct.atimp1[i-1] = 0;
         pstruct.atimp2[i-1] = 0;
         pstruct.atimp3[i-1] = 0;
         pstruct.atimp4[i-1] = 0;
      }
   }
 
   for(i=FirstNBond; i<=values.nnbs; i++)
   {
      InSegment = (pstruct.nbexcl[i-1] >= FirstAtom && 
                   pstruct.nbexcl[i-1] <= values.natoms) ?
                  TRUE : FALSE;
      if(!InSegment) pstruct.nbexcl[i-1] = 0;
   }
}
 

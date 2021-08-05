#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
void updtnictot(
int mode,
int prflag
)
{
   makinb(mode);
 
   if(prflag) fprintf(out,"Current internal coordinate counts:\n\n\
   NATOM   NBOND  NTHETA    NPHI  NIMPHI     NNB\n%8d%8d%8d%8d%8d%8d\n",
values.natoms,values.nbonds,values.nangs,values.nptors,values.nitors,
values.nnbs);
 
   pstruct.segndx[values.nsegs][0] = values.nres;
   pstruct.segndx[values.nsegs][1] = values.natoms;
   pstruct.segndx[values.nsegs][2] = values.nbonds;
   pstruct.segndx[values.nsegs][3] = values.nangs;
   pstruct.segndx[values.nsegs][4] = values.nptors;
   pstruct.segndx[values.nsegs][5] = values.nitors;
   pstruct.segndx[values.nsegs][6] = values.nnbs;
   pstruct.segndx[values.nsegs][7] = values.ndonat;
   pstruct.segndx[values.nsegs][8] = values.naccat;
}
 

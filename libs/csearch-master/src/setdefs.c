#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Setup various default values */
 
void setdefs(
void
)
{
   /* Set all the FORTRAN unit numbers to -1 to catch bugs */
   cg.glymapu  = -1;
   cg.alamapu  = -1;
   cg.promapu  = -1;
   cg.proconsu = -1;
   cg.stunit   = -1;
   /* Now the energies */
   cg.glyemax  = 100.0;
   cg.alaemax  = 100.0;
   cg.proemax  = 100.0;
   cg.eringpro = largnum;
   /* Other bits and pieces */
   cg.ignore_evdw = f77_false;
   dbg.cgen = 0;
   cg.maxleaf = largint;
   cg.restart_st = NULL;
   cg.restart_stlen = 0;
   cg.save_coor = f77_true;
   cg.sidehits_opt = f77_true;
   /* Non-bonded stuff */
   engpar.nbcut = 5.0;
   engpar.dielec = 50.0;
   hbonds.hbcut= 4.5;
   hbonds.hbacut = 90.0;
   /* RTF max array sizes */
   restop.nparmx[0] = mxatmr;    /* Atoms       */
   restop.nparmx[1] = mxbndr;    /* Bonds       */
   restop.nparmx[2] = mxangr;    /* Angles      */
   restop.nparmx[3] = mxtorr;    /* Torsions    */
   restop.nparmx[4] = mximpr;    /* Impropers   */
   restop.nparmx[5] = mxnber;    /* Exclusions  */
   restop.nparmx[6] = mxhbdr;    /* Donors      */
   restop.nparmx[7] = mxhbar;    /* Acceptors   */
   restop.nparmx[8] = mxicr;     /* Builds      */
   restop.nparmx[9] = mxgrpr;    /* Groups      */
   /* Zero the current counts */
   values.nres   = 0;
   values.natoms = 0;
   values.nbonds = 0;
   values.nangs  = 0;
   values.nptors = 0;
   values.nitors = 0;
   values.nnbs   = 0;
   values.nsegs   = 0;
   /* The first position in nictot is always 0 */
   pstruct.segndx[0][0] = 0;
   pstruct.segndx[0][1] = 0;
   pstruct.segndx[0][2] = 0;
   pstruct.segndx[0][3] = 0;
   pstruct.segndx[0][4] = 0;
   pstruct.segndx[0][5] = 0;
   pstruct.segndx[0][6] = 0;
   pstruct.segndx[0][7] = 0;
   pstruct.segndx[0][8] = 0;
   pstruct.segndx[0][9] = 0;
}

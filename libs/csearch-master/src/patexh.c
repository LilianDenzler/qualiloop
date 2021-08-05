#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/* Patch the internal coordinates for coordinates specific to proteins
   using explicit hydrogens.
   Assumes all residues are amino acids. Will screw up, if they're not.

   02.08.92 Some rewriting.   By:   ACRM.
*/
 
#define ATOM_C    "C   "
#define ATOM_CA   "CA  "
#define ATOM_OT1  "OT1 "
#define RES_PRO   "PRO "
#define ATOM_HC   "HC  "
#define ATOM_NH3  "NH3 "
#define ATOM_OC   "OC  "
#define ATOM_NT   "NT  "
#define ATOM_HT3  "HT3 "
#define ATOM_ACET "CH3 "
 
void patexh(
int ResNum,
int AtomNum,
int TorsionNum,
int ImproperNum,
int DonorNum
)
{
   int   i, 
         FirstC,
         HC_ptr,  NH3_ptr, OC_ptr,
         i_value, 
         TorsionCount, AtomCount;
   float f_value;

   /* Find indexes of HC, NH3 and OC in the residue topology information */
   for(i=1; i<=values.natyps; i++)
   {
      if(!strncmp(restop.acodes[i-1],ATOM_HC, 4)) HC_ptr  = i;
      if(!strncmp(restop.acodes[i-1],ATOM_NH3,4)) NH3_ptr = i;
      if(!strncmp(restop.acodes[i-1],ATOM_OC, 4)) OC_ptr  = i;
   }
   
   /* Find the first C atom */
   FirstC = matom(ResNum+1,ATOM_C);

   /* Providing it's not residue PCA, and the atom is not CH3, we need
      to patch the N-terminus
   */
   if(strncmp(pstruct.atmnme[AtomNum-1],ATOM_ACET,4) &&
      strncmp(pstruct.resnme[ResNum-1],"PCA ",4))
   {
      /* Modify Nter internal coords */
      strncpy(pstruct.atmnme[AtomNum+1],ATOM_NT,4);
      pstruct.atimp1[ImproperNum-1] = 0;
      if(!strncmp(pstruct.resnme[ResNum],RES_PRO,4))  /* Proline                 */
      {
         pstruct.attor3[TorsionNum-1]--;
         pstruct.attor4[TorsionNum-1] = matom(ResNum+1,ATOM_C);
         pstruct.attor1[TorsionNum]   = 0;
         pstruct.atcode[AtomNum+1]   = NH3_ptr;
      }
      else                                      /* Not proline             */
      {
         pstruct.attor4[TorsionNum-1] = matom(ResNum+1,ATOM_C);
         pstruct.attor4[TorsionNum]   = pstruct.attor4[TorsionNum-1];
         pstruct.atcode[AtomNum+1]   = NH3_ptr;
         pstruct.atcode[AtomNum+2]   = HC_ptr;
         strncpy(pstruct.atmnme[AtomNum+2],ATOM_HT3,4);
         f_value = (pstruct.atchrg[AtomNum-1] + pstruct.atchrg[AtomNum] 
                  + pstruct.atchrg[AtomNum+2])/3.0;
         pstruct.atchrg[AtomNum-1]   = pstruct.atchrg[AtomNum]
                                     = pstruct.atchrg[AtomNum+2] = f_value;
      }
      
      i_value = matom(ResNum+1,ATOM_CA);
      pstruct.atchrg[i_value-1] += 0.020;      /* Patch charge on Nter CA    */
      pstruct.atchrg[AtomNum+1] += 0.261;      /* Patch charge on Nter N     */

      for(i=1; i<=3; i++)
      {
         /* Fix the HBond donors, since the N-ter will be NH3 not NH */
         if(strncmp(pstruct.resnme[ResNum],RES_PRO,4) || i!=3)
         {
            if(!strncmp(pstruct.resnme[ResNum],RES_PRO,4))   /* Proline */
               pstruct.hbdan1[DonorNum+i-2] = AtomNum + 3;
            else                                             /* Not proline  */
               pstruct.hbdan1[DonorNum+i-2] = AtomNum + 4;
            pstruct.hbdan2[DonorNum+i-2] = FirstC;
            pstruct.hbdonr[DonorNum+i-2] = AtomNum + 2;
         }
      }
 
      /* Modify Cter internal coordinates */
      i_value = values.nptors + 3 - 
                restop.nparam[pstruct.resndx[values.nres-2]-1][3];
      pstruct.attor1[i_value-1]       = 0;
      pstruct.nbexcl[values.nnbs-4]   = 0;
      pstruct.nbexcl[values.nnbs-3]   = 0;
      pstruct.atcode[values.natoms-2] = OC_ptr;
      pstruct.atchrg[values.natoms-2] -= 0.200;   /* Patch charge on Cter OT1 */
      pstruct.atchrg[values.natoms-3] += 0.030;   /* Patch charge on Cter C   */
      strncpy(pstruct.atmnme[values.natoms-2],ATOM_OT1,4);
 
      /* Mods for prolines. TorsionCount keeps track of previous residue   */
      TorsionCount = TorsionNum-1;
      for(i=ResNum+1; i<=values.nres; i++)
      {
         if(!strncmp(pstruct.resnme[i-1],RES_PRO,4))
         {
            AtomCount = pstruct.lstatm[i-1];
            /* Check if first res in a segment is a PRO as exclusions are
               different for the Nter. Also a torsion coming into the proline
               must be fixed.
            */
            if(i==ResNum+1)
            {
               pstruct.nbexcl[pstruct.lstexc[AtomCount-1]-1] += 2;
               pstruct.nbexcl[pstruct.lstexc[AtomCount-1]-4] += 2;
            }
            else
            {
               pstruct.nbexcl[pstruct.lstexc[AtomCount-1]-2] += 2;
               pstruct.attor4[TorsionCount+2]--;
            }
         }
         TorsionCount += restop.nparam[pstruct.resndx[i-2]-1][3];
      }
   }
}
 

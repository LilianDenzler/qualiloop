#include "ProtoTypes.h"
#include "CongenProto.h"
 
void fillbbatms(
struct atom **atoms,
int resno,
logical forward,
logical nter,
struct atom **at_n,
struct atom **at_h,
struct atom **at_ca,
struct atom **at_cb,
struct atom **at_cg,
struct atom **at_cd,
struct atom **at_c,
struct atom **at_o,
struct atom **at_ht1
)
{
   int firstat,lastat;
 
   firstat = pstruct.lstatm[resno-1] + 1;
   lastat = pstruct.lstatm[resno];
 
   if (forward)
   {
      fillatm_ca(at_ca,resno,forward);
      fillatm_c(at_c,resno,forward);
      fillatm_o(at_o,resno,forward);
      fillatm_n(at_n,resno+1,forward);
   }
   else
   {
      fillatm_n(at_n,resno,forward);
      fillatm_c(at_c,resno-1,forward);
      fillatm_o(at_o,resno-1,forward);
      fillatm_ca(at_ca,resno-1,forward);
   }
   *at_cb = null;
   *at_cg = null;
   *at_cd = null;
   *at_h = null;
   *at_ht1 = null;
   if (strncmp(pstruct.resnme[resno-1],"PRO ",4) == 0)
   {
      fillatm_cb(at_cb,resno,true);
      *at_cg = (struct atom *) alloc(sizeof(struct atom));
      (*at_cg)->atomno = srwdbd(pstruct.atmnme,&firstat,&lastat,"CG  ");
      clear_atom(*at_cg);
      *at_cd = (struct atom *) alloc(sizeof(struct atom));
      (*at_cd)->atomno = srwdbd(pstruct.atmnme,&firstat,&lastat,"CD  ");
      clear_atom(*at_cd);
   }
   else if (strncmp(pstruct.resnme[resno-1],"GLY ",4) == 0)
   {
      fillatm_h(at_h,resno,forward);
   }
   else
   {
      fillatm_h(at_h,resno,forward);
      fillatm_cb(at_cb,resno,true);
   }
 
   if(nter && !forward) fillatm_ht1(at_ht1,resno-1,false);
   if (*at_n != null) *atoms++ = *at_n;
   if (*at_ca != null) *atoms++ = *at_ca;
   if (*at_c != null) *atoms++ = *at_c;
   if (*at_o != null) *atoms++ = *at_o;
   if (*at_cb != null) *atoms++ = *at_cb;
   if (*at_cg != null) *atoms++ = *at_cg;
   if (*at_cd != null) *atoms++ = *at_cd;
   if (*at_h != null) *atoms++ = *at_h;
   if (*at_ht1 != null) *atoms++ = *at_ht1;
   *atoms = null;
}
 

#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Check for any bad contacts to the atoms in the list.
*/
 
logical chkcntcts(
struct atom **atomsp,
float *maxevdw,
int *sidehits
)
{
   logical ok;
   logical impact;
   struct atom **p,**q;
   int mode;
   double eel,evdw,ehb;
   char buf[21];
 
   mode = contact;
   p = atomsp;
   ok = true;
   for (p = atomsp; *p; p++)
   {
      /* ACRM added floats so schnatm is called properly */
      float eel_f  = (float)eel,
            evdw_f = (float)evdw,
            ehb_f  = (float)ehb;
 
      schnatm((*p)->atomno,mode,maxevdw,cg.ignore_evdw,sidehits,
                       &impact,&eel_f,&evdw_f,&ehb_f);
      eel  = (double)eel_f;
      evdw = (double)evdw_f;
      ehb  = (double)ehb_f;
 
      if (impact)
      {
         ok = false;
         if (dbg.cgen)
         {
            buf[0] = '\0';
            print_atom1(buf,(*p)->atomno);
            fprintf(out,"Check_contacts: Bad contact for atom %s\n",buf);
         }
         break;
      }
      adatmtgrd((*p)->atomno);
   }
   if (!ok)
      for (q = atomsp; q < p; q++) delatmfgrd((*q)->atomno);
   return ok;
}
 

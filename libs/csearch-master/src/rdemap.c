#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Reads in the energy map for a particular peptide. */
 
void rdemap(
int unit,
int *nmap,
struct pepmap **map,
float emax,
char name[]
)
{
   FILE *cunit;
   float emin,emincis,e1,e2,omega,phi,psi;
   int i,nentries,count;
#define LINEMX 200
   int linemx = LINEMX;
   char line[LINEMX];
#undef LINEMX
   struct pepmap *p;
   logical ok,cis_ok;
 
   if(unit==1)
      cunit = glymap_fp;
   else if(unit==2)
      cunit = promap_fp;
   else if(unit==3)
      cunit = alamap_fp;
   else
   {
      fprintf(out,"Unknown peptide map unit in rdemap: %d\n",unit);
      die();
   }
 
   rewind(cunit);
/* We don't do titles properly
//   rdtitl(ctitla.titlea,&ctitla.ntitla,&unit);
*/
 
   fprintf(out,"\nReading energy map for %s\n",name);
   rdprtitle(cunit);
 
/*
//   prtitl(ctitla.titlea,&ctitla.ntitla,(i=6,&i));
*/
   if (!fgets(line,linemx,cunit))
   {
      fprintf(out,"\nError in RDEMAP for map %s\nUnexpected EOF 1\n",name);
      die();
   }
   terminate(line);
   sscanf(line,"%d",&nentries);
   if (nentries <= 0) nentries = largint;
   emin = largnum;
   count = 0;
   emincis = largnum;
   while (count < nentries)
   {
      if (!fgets(line,linemx,cunit)) break;
      terminate(line);
      count++;
      scnemapln(line,&omega,&phi,&psi,&e1,&e2);
      emin = minf2(emin,e1);
      emincis = minf3(emincis,e1,e2);
   }
   fprintf(out,"Minimum energy is %g; minimum cis peptide energy is %g\n",
          emin,emincis);
   rewind(cunit);
/*
//   rdtitl(ctitla.titlea,&ctitla.ntitla,&unit);
*/
   skiptitle(cunit);
   if (!fgets(line,linemx,cunit))
   {
      fprintf(out,"\nError in RDEMAP for map %s\nUnexpected EOF 2\n",name);
      die();
   }
   terminate(line);
   *map = (struct pepmap *) alloc(count*sizeof(struct pepmap));
   for (p = *map, i=0; i < count; i++)
   {
      if (!fgets(line,linemx,cunit))
      {
         fprintf(out,"\nError in RDEMAP for map %s\nUnexpected EOF 2\n",name);
         die();
      }
      terminate(line);
      scnemapln(line,&omega,&phi,&psi,&e1,&e2);
      ok  = e1-emin < emax && fabs(omega) >= 90.0;
      cis_ok = e1-emin <= emax || e2-minf2(emin,emincis) <= emax;
      if (ok || cis_ok)
      {
         p->omega = omega * dtorad;
         p->phi = phi * dtorad;
         p->psi = psi * dtorad;
         if (ok)
         {
            p->energy = e1;
            p->cisonly = false;
         }
         else
         {
            p->energy = minf2(e1,e2);
            p->cisonly = true;
         }
         p++;
      }
   }
   *nmap = p - *map;
   fprintf(out,"%d entries were selected.\n",*nmap);
}
 

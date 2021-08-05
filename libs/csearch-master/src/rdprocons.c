#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Reads in the proline constructor arrays. */
 
void rdprocons(
int unit,
float (**procons)[9],
float **proconsphi,
float (**eprocons)[3],
int *nprocons,
float eringpro
)
{
   int i,j,n;
#define LINEMX 200
   char line[LINEMX];
   int linemx = LINEMX;
#undef LINEMX
   float emin,e,ec,*pc,*pcp,*ep,pca[9],pcpa,epa[3];
 
   rewind(procons_fp);
/* We don't handle titles
//   rdtitl(ctitla.titlea,&ctitla.ntitla,&unit);
*/
   fprintf(out,"\nTitle for proline construction file:\n\n");
   rdprtitle(procons_fp);
/*
//   prtitl(ctitla.titlea,&ctitla.ntitla,(i=6,&i));
*/
   if (!fgets(line,linemx,procons_fp))
   {
      fprintf(out,"\nError in RDPROCONS -- Unexpected EOF 1\n");
      die();
   }
   terminate(line);
   sscanf(line,"%d",&n);
   emin = largnum;
   for (i=0; i<n; i++)
   {
      if (!fgets(line,linemx,procons_fp))
      {
          fprintf(out,"\nError in RDPROCONS -- Unexpected EOF 2\n");
          die();
      }
      terminate(line);
      sscanf(line,"%*10c %10f %*10c %10f",&e,&ec);
      emin = minf2(emin,e-ec);
   }
   fprintf(out,"\nMinimum ring energy is %g\n",emin);
   rewind(procons_fp);
/*
//   rdtitl(ctitla.titlea,&ctitla.ntitla,&unit);
*/
   skiptitle(procons_fp);
 
   if (!fgets(line,linemx,procons_fp))
   {
      fprintf(out,"\nError in RDPROCONS -- Unexpected EOF 3\n");
      die();
   }
   terminate(line);
   *procons = (float *) alloc(n*sizeof(float)*9);  /* ACRM 20.02.06 FIXME */
   *proconsphi = (float *) alloc(n*sizeof(float)); /* ACRM 20.02.06 FIXME */
   *eprocons = (float *) alloc(n*sizeof(float)*3); /* ACRM 20.02.06 FIXME */
   for (pcp = *proconsphi, pc = *procons, ep = *eprocons, *nprocons = 0, i = 0;
        i<n;
        i++)
   {
      if (!fgets(line,linemx,procons_fp))
      {
         fprintf(out,"\nError in RDPROCONS -- Unexpected EOF 4\n");
         die();
      }
      terminate(line);
      sscanf(line,"%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f",
             &pcpa,&epa[0],&epa[1],&epa[2],
             &pca[0],&pca[1],&pca[2],&pca[3],&pca[4],
             &pca[5],&pca[6],&pca[7],&pca[8]);
      if (epa[0]-epa[2] <= emin+eringpro)
      {
         *pcp++ = pcpa*dtorad;
         for (j=0; j<3; j++) *ep++ = epa[j];
         *pc++ = pca[0];
         *pc++ = pca[1]*dtorad;
         *pc++ = pca[2]*dtorad;
         *pc++ = pca[3];
         *pc++ = pca[4]*dtorad;
         *pc++ = pca[5]*dtorad;
         *pc++ = pca[6];
         *pc++ = pca[7]*dtorad;
         *pc++ = pca[8]*dtorad;
         (*nprocons)++;
      }
   }
   fprintf(out,"\n%d ring constructors selected\n",*nprocons);
   fprintf(out,"allowing phi angles between %g and %g\n",
          **proconsphi/dtorad,*--pcp/dtorad);
}
 

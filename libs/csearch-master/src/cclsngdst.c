#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*
*   Given a list of atoms in catom (natom in number) for the backbone
*   degree of freedom pointed to by desc, finds the maximum possible
*   closing distance. If stretch_all is true, then all bond angles are
*   stretched by cmaxdt[0]. Otherwise, any angle appearing in cangle is
*   adjusted by the corresponding entry in cmaxdt.
*
*   A list of four coordinates is kept with the first atom never changed,
*   but the last three constantly being shuffled. CARTX2 is used for
*   successive constructions, the torsion angle always being trans.*/
 
 
void cclsngdst (
struct backbone_d *desc,
int *catom,
int natom,
logical stretch_all,
int **cangle,
float cmaxdt[],
int nangle
)
{
   float x[4],y[4],z[4];
   float b,t;
   int i,ires;
   char buf[21];
 
/* COM OMUPD BNJ */
/* local version of pi to be sent to cartx2...!!! */
   float phi = PI;
 
   if (dbg.cgen > 0)
   {
      ires = desc->resno-1;
      fprintf(out,"Atoms for closing distance for residue %.4s %.4s\n",
             pstruct.resnme[ires],pstruct.resid[ires]);
      for (i=0; i<natom; i++)
      {
         buf[0] = '\0';
         print_atom1(buf,catom[i]);
         fprintf(out,"%s\n",buf);
      }
   }
   if (natom <= 1)
   {
      desc->closing_distance = 0.0;
      return;
   }
   for (i=0; i<4; i++) z[i] = 0.0;
   b = getparbond2(catom[0],catom[1]);
   x[0] = 0.0;
   y[0] = 0.0;
   x[1] = b;
   y[1] = 0.0;
   if (natom == 2)
   {
      desc->closing_distance = b;
      return;
   }
   b = getparbond2(catom[1],catom[2]);
   t = gtprangclsa(catom[0],catom[1],catom[2],stretch_all,
                        cangle,cmaxdt,nangle);
   x[2] = x[1] - b*(float)cos((double)t);
   y[2] = b*(float)sin((double)t);
   if (natom == 3)
   {
      desc->closing_distance = sqrt(x[2]*x[2] + y[2]*y[2]);
      return;
   }
   for (i=3;;)
   {
      b = getparbond2(catom[i],catom[i-1]);
      t = gtprangclsa(catom[i],catom[i-1],catom[i-2],
                           stretch_all,cangle,cmaxdt,nangle);
 
      cartx2(x,y,z,&one,&two,&three,&four,&b,&t,&phi);
 
      if (++i >= natom) break;
      x[1] = x[2];
      y[1] = y[2];
      z[1] = z[2];
      x[2] = x[3];
      y[2] = y[3];
      z[2] = z[3];
   }
   desc->closing_distance = sqrt(x[3]*x[3] + y[3]*y[3]);
   if (dbg.cgen > 0)
      fprintf(out,"Z error in CALCULATE_CLOSING_DISTANCE is %g\n",z[3]);
   return;
}
 

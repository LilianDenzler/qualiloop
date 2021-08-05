#include "ProtoTypes.h"
#include "CongenProto.h"
 
#ifndef PI
#define PI ((double)4.0*atan((double)1.0))
#endif
#define TWOPI ((double)2.0 * PI)
#define RADI  ((double)180.0 / PI)
#define ANUM 9999.0

#ifdef TORSION
#undef TORSION
#endif
#define TORSION   1

#define IMPROPER  2

#define TAYLORCUT 0.1

extern int wrnlev;
extern FILE *out;
 
/****************************************************************************
 Routine  - ephi()

 Function - This routine returns the potential energy due to proper and
            improper torsion angles. The energy is expressed as a function
            of the cosine of the dihedral angle rather than the angle itself
            to avoid continuity problems, particularly when the dihedrals
            become planar.

 Equations                                                             
  proper  - Energy = Constant * { 1 + phase * F(periodicity,cosine) }
            where:
                  cosine = cosine(dihedral angle)
                  phase  = +1 or -1
                  F(1,cosine) = cosine
                  F(2,cosine) = 2*(cosine**2)-1
                  F(3,cosine) = 4*(cosine**3)-3*cosine 
                  F(4,cosine) = 8*(cosine**4)-8*(cosine**2)+1
                  F(5,cosine) = error
                  F(6,cosine) = 32*(cosine**6)-48*(cosine**4)+18*(cosine**2)-1 

 improper - Energy = Constant * (phi - phiref)**2
            where:
                  phi   = dihedral angle
                  phiref= a reference angle

            If phiref = 0 or abs(sin(phi) < TAYLORCUT the equation is 
            expanded into a 3rd order Taylor series in cosine(phi).
 
   31.07.92 Rewritten   Author:  ACRM
   09/11/92 rkw changed the name of Angle to AngVal
****************************************************************************/
 
double   ephi(int       Type,          /* 1: Torsions, 2: Impropers        */
              short int *AtArray1,     /* Arrays of indexes to the atoms   */
              short int *AtArray2,  
              short int *AtArray3,
              short int *AtArray4,
              short int *TorsionIdx,   /* Array of indexes to the torsions */
              int       NumPhis,       /* Number of torsions               */
              float     *ForceCon,     /* Arrays of parameters for phis    */
              float     *Period,
              float     *OptAngle,
              float     *x,            /* Coordinate arrays                */
              float     *y,
              float     *z,
              short int *Fixed,        /* Flags indicating immobile atoms  */
              int       PrintFlag,     /* Flag to print info               */
              int       MaxPrint)      /* Max number of items to print     */
{
   double   TotEnergy;
   float    Energy,     Force,      DeltaAngle,
            SinAngle,   CosAngle,   CosAngleSq,  
            AngVal,
            arg,        AngDegrees, OptDegrees;
   int      NBentWarnings,
            NLinearWarnings,
            PhiCount,
            TorNumber,
            i,    j,    k,    l,
            IntPeriod;

   
   if(PrintFlag) fprintf(out,"\n\n   PHI  ATOM ATOM ATOM ATOM  ENERGY\
    ANGLE     FORCE-CON DELTA MINIMUM\n");
    
   TotEnergy      = 0.0;
   if(NumPhis==0) return(TotEnergy);

   NBentWarnings   = 0;
   NLinearWarnings = 0;
   
   for(PhiCount=0; PhiCount<NumPhis; PhiCount++)
   {
      i = AtArray1[PhiCount] - 1;
      if(i<0) continue;

      j = AtArray2[PhiCount] - 1;
      k = AtArray3[PhiCount] - 1;
      l = AtArray4[PhiCount] - 1;
      
      /* Skip this torsion if any atom == ANUM, or all atoms imobile */
      if(x[i]==ANUM || 
         x[j]==ANUM || 
         x[k]==ANUM || 
         x[l]==ANUM ||
         (Fixed[i] && Fixed[j] && Fixed[k] && Fixed[l]))
         continue;
         
      TorNumber = TorsionIdx[PhiCount] - 1;
      
      /* Check and warn for linearity */
      if((fabs((double)atomangle(x[i],y[i],z[i],x[j],y[j],z[j],
                                 x[k],y[k],z[k])) <= 0.1) ||
         (fabs((double)atomangle(x[j],y[j],z[j],x[k],y[k],z[k],
                                 x[l],y[l],z[l])) <= 0.1))
      {
         if((NLinearWarnings++<5 && wrnlev>=0) || wrnlev>0)
            fprintf(out,"Warning==> Phi %d is almost linear while calculating \
torsion energy.\nAtoms: %d %d %d %d\n",TorNumber+1,i+1,j+1,k+1,l+1);
      }

      /* Calculate the torsion angle */
      AngVal = phi(x[i],y[i],z[i],x[j],y[j],z[j],
                  x[k],y[k],z[k],x[l],y[l],z[l]);

      SinAngle  = (float)sin((double)AngVal);
      CosAngle  = (float)cos((double)AngVal);

/* Energy contribution */

      if(Type==TORSION)         /* Proper dihedrals */
      {

         IntPeriod = (int)(Period[TorNumber]+0.0001);
         
         if(IntPeriod==(int)(Period[TorNumber]-0.0001))  goto error;
         if(IntPeriod>6 || IntPeriod<1)                  goto error;
         
         switch(IntPeriod)
         {
         case 1:
            Energy     = CosAngle;
            break;
         case 2:
            Energy     = 2.0*CosAngle*CosAngle - 1.0;
            break;
         case 3:
            CosAngleSq = CosAngle * CosAngle;
            Energy     = CosAngle * (4.0*CosAngleSq - 3.0);
            break;
         case 4:
            CosAngleSq = CosAngle * CosAngle;
            Energy     = 1.0 + CosAngleSq*8.0*(CosAngleSq-1.0);
            break;
         case 5:
            goto error;
         case 6:
            CosAngleSq = CosAngle  * CosAngle;
            Energy     = CosAngleSq * 
                         (CosAngleSq*(CosAngleSq*32.0 - 48.0) + 18.0) - 1.0;
            break;
         }
         
         arg = ForceCon[TorNumber];
         if(OptAngle[TorNumber]!=0.0)
         {
            arg = (-arg);
            if(fabs(PI-OptAngle[TorNumber])>0.01) goto error;
         }
         Energy     = ForceCon[TorNumber] + arg*Energy;
         
         TotEnergy += (double)Energy;
         
         if(PrintFlag && PhiCount<MaxPrint)
         {
            AngDegrees = RADI * AngVal;
            OptDegrees = RADI * OptAngle[TorNumber];
            fprintf(out,"%5d%5d%5d%5d%5d%10.4f%10.4f%12.4f%6.2f%9.4f\n",
                    PhiCount+1,i+1,j+1,k+1,l+1,Energy,AngDegrees,
                    ForceCon[TorNumber],Period[TorNumber],OptDegrees);
         }
      }
      else  /* Type == IMPROPER, Improper dihedrals */
      {
         DeltaAngle = AngVal-OptAngle[TorNumber];
         if(DeltaAngle >  PI) DeltaAngle -= TWOPI;
         if(DeltaAngle < -PI) DeltaAngle += TWOPI;

         if(OptAngle[TorNumber]==0.0)
         {
            /* TAYLORCUT can be modified to minimise error */
            if(fabs(SinAngle) > TAYLORCUT)
            {
               Force   = DeltaAngle * ForceCon[TorNumber];
               Energy  = DeltaAngle * Force;
            }
            else
            {
               Energy  = AngVal * AngVal * ForceCon[TorNumber];
            }
         }
         else
         {
            AngDegrees = RADI * AngVal;
            if(fabs(SinAngle) < 0.001) SinAngle=0.001;
            if(fabs(DeltaAngle) >= PI/2.0)
            {
               NBentWarnings++;
               if((NBentWarnings<=5 && wrnlev>=0) || wrnlev>=1)
                  fprintf(out,"Warning==> Bent improper torsion angle \
is too far\nfrom mimimum for IPHI=%d and PHI=%14.8f while calculating \
energy.\nAtoms %5d %5d %5d %5d\n",
PhiCount+1,AngDegrees,i+1,j+1,k+1,l+1);
            }
            Force   = DeltaAngle * ForceCon[TorNumber];
            Energy  = DeltaAngle * Force;
         }
         TotEnergy += (double)Energy;

         if(PrintFlag && PhiCount<MaxPrint)
         {
            AngDegrees = RADI * AngVal;
            OptDegrees = RADI * OptAngle[TorNumber];
            fprintf(out,"%5d%5d%5d%5d%5d%10.6f%10.4f%12.4f%15.4f\n",
                    PhiCount+1,i+1,j+1,k+1,l+1,Energy,AngDegrees,
                    ForceCon[TorNumber],OptDegrees);
         }
      }
   }
      
      
   if((NBentWarnings += NLinearWarnings) > 5)
      fprintf(out,"Total of %6d warnings while calculating torsion energy\n",
              NBentWarnings);

   return(TotEnergy);

error:
   fprintf(out,"\nError(s) in parameter list for dihedral angles \
while calculating energy.\n\n");
   die();
   
   return(TotEnergy);
}


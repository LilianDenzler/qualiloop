#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Does a Quicksort on an array
   Recoded into C from FORTRAN in Numerical Recipes */
 
#define QSORT_CUT 7
#define RAND_FM 7875.0
#define RAND_FA 211.0
#define RAND_FC 1663.0
#define RAND_FMI 1.0/RAND_FM
 
void sort_ui(
unsigned long int *arr,
int n
)
{
   int *istack;
   int jstack = 0;
   int l = 1;
   int ir = n;
   float fx = 0.0;
   int i, j, iq,
       nstack;
 
   unsigned long int a;
 
   /* Work out the required stack space and allocate it */
   nstack = 2 * (int)(log((double)n) / log((double)2.0)) + 1;
   istack = (int *)malloc(nstack * sizeof(int));
 
   for(;;)
   {
      /* If it's small just do a straight insertion sort */
      if(ir-l < QSORT_CUT)
      {
         for(j= l+1; j<=ir; j++)
         {
            a = arr[j-1];
            for(i=j-1; i>=1; i--)
            {
               if(arr[i-1] <= a) break;
               arr[i] = arr[i-1];
            }
            arr[i] = a;
         }
         /* Free the stack and exit when the stack is empty */
         if(jstack==0)
         {
            free(istack);
            return;
         }
 
         ir = istack[jstack-1];
         l = istack[jstack-2];
         jstack -= 2;
      }
      else /* Do the partition sort */
      {
         i=l;
         j = ir;
 
         /* Randomise the starting order */
         fx = (float)fmod((double)(fx*RAND_FA+RAND_FC),(double)RAND_FM);
         iq = l + (ir-l+1)*(fx*RAND_FMI);
 
         a = arr[iq-1];
         arr[iq-1] = arr[l-1];
         for(;;)
         {
            while(j>0 && a<arr[j-1]) j--;
            if(j<=i)
            {
               arr[i-1] = a;
               break;
            }
            arr[i-1] = arr[j-1];
            i++;
            while(i<=n && a>arr[i-1]) i++;
            if(j<=i)
            {
               arr[j-1] = a;
               i = j;
               break;
            }
            arr[j-1] = arr[i-1];
            j--;
         }
         jstack += 2;
         if(jstack > nstack)
         {
            printf("sort_i(): Allocation of nstack not big enough!\n");
            exit(1);
         }
         if(ir-i >= i-l)
         {
            istack[jstack-1] = ir;
            istack[jstack-2] = i+1;
            ir = i-1;
         }
         else
         {
            istack[jstack-1] = i-1;
            istack[jstack-2] = l;
            l = i+1;
         }
      }
   }
   return;
}
 

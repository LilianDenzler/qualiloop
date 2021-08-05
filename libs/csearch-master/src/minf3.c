#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Computes minimum of three floating point numbers. */
 
float minf3(
float x1,
float x2,
float x3
)
{
   float m;
 
   m = x1 < x2 ? x1 : x2;
   m = m < x3 ? m : x3;
   return m;
}
 

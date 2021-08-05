/* Generally useful macros */

#ifndef USHORT
#define USHORT unsigned short
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define NEXT(z) (z)=(z)->next
#define ALLOCNEXT(z,y) (z)->next = (y *)malloc(sizeof(y)); (z)=(z)->next
#define FREELIST(y,z) while((y)) \
                       { z *q; \
                         q = (y)->next; \
                         free((char *)(y)); \
                         (y) = q; \
                       }
                          
#define CR 13
#define LF 10

typedef struct
{
   float x,y,z;
}  Vec3f;

#define DISTSQ(a,b) (((a)->x - (b)->x) * ((a)->x - (b)->x) + \
                     ((a)->y - (b)->y) * ((a)->y - (b)->y) + \
                     ((a)->z - (b)->z) * ((a)->z - (b)->z))

#define DIST(a,b) sqrt(((a)->x - (b)->x) * ((a)->x - (b)->x) + \
                       ((a)->y - (b)->y) * ((a)->y - (b)->y) + \
                       ((a)->z - (b)->z) * ((a)->z - (b)->z))

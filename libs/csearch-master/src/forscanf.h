#ifndef __FORSCANF_H__
#define __FORSCANF_H__

struct flist
{
    struct flist *next;
    int   type;
    int   retd;
    float retf;
    char  *rets;
    char  retc;
}   ;
typedef struct flist FLIST;
 
#define FREEFORSCANFLIST(y)   while((y)) \
                              {  FLIST *q; \
                                 if((y)->type == 3 && (y)->rets) free((y)->rets); \
                                 q = (y)->next; \
                                 free(y); \
                                 (y) = q; \
                              }

#endif

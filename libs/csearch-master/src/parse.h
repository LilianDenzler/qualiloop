#ifndef __PARSE_H__
#define __PARSE_H__

 
/* Defines used for the MAKEKEY macro */
#define STRING 1
#define NUMBER 0
 
/* Return values from parse() */
#define PARSE_ERRC    -1
#define PARSE_ERRP    -2
#define PARSE_COMMENT -3
 
/* Definition of the KeyWd structure in which we store keywords */
typedef struct
{
   char  *name;
   int   string, nparam;
}  KeyWd;
 
/* Macro to create a keyword */
#define MAKEKEY(x,w,v,z) \
        (x).name = (char *)malloc((strlen(w)+2) * sizeof(char)); \
        strcpy((x).name,w); \
        (x).string = v; \
        (x).nparam = z

#endif

#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Searches the words array for 4-letter word */
 
int srchwd(
char words[][4],
int  numwrd,
char word[]
)
{
int i;
 
   if(numwrd == 0) return(0);
   for(i=0;i<numwrd;i++)
      if(!strncmp(word,words[i],4))
          return(i+1);
   return(0);
}
 

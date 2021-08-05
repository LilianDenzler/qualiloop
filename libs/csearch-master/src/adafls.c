#include "ProtoTypes.h"
#include "CongenProto.h"

/* -- FUNCTION -- */
/* This is used for linked list handling. It is called from FORTRAN, so
   uses pointers in the function call.
   
   02.08.92 Recoded.    By:   ACRM

*/
int adafls(int *ListStart,
           int *ListEnd,
           int *Position,
           int *NewItem,
           int *NextItem)
{
   if(*Position == 0)
   {
      if(*ListStart == 0)                 /* Initialise list            */
      {
         *ListStart = *NewItem;
         *ListEnd   = *NewItem;
         NextItem[*NewItem-1] = 0;
      }
      else                                /* Insert at start of list    */
      {
         NextItem[*NewItem-1] = *ListStart;
         *ListStart = *NewItem;
      }
   }
   else if(NextItem[*Position-1] != 0)    /* Middle of list             */
   {
      NextItem[*NewItem-1]  = NextItem[*Position-1];
      NextItem[*Position-1] = *NewItem;
   }
   else if(*ListEnd!=*Position)           /* List or position is fubar  */
   {
      fprintf(out,"Bad position passed to adafls list handler.\n");
      die();
   }
   else                                   /* End of list                */
   {
      NextItem[*ListEnd-1] = *NewItem;
      NextItem[*NewItem-1] = 0;
      *ListEnd             = *NewItem;
   }
   
   return(0);
}

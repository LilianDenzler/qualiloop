/*************************************************************************

   Program:    
   File:       fsscanf.c
   
   Version:    V1.1
   Date:       12.07.93
   Function:   Read from a string using FORTRAN-like rigid formatting
   
   Copyright:  (c) SciTech Software 1993
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0372) 275775
   EMail:      UUCP:  {...}cbmehq!cbmuk!scitec!amartin
                      amartin@scitec.adsp.sub.org
               JANET: andrew@uk.ac.ox.biop
               
**************************************************************************

   This code is copyright (c) SciTech Software, 1993. It is not in the 
   public domain. In the future it is hoped that it will form part of a 
   book on manipulating protein structures in C; along the lines of
   `Numerical Recipes'.
   
   It may be freely copied and distributed ONLY as part of the QTree
   CPK picture program (and associated programs) from SciTech Software. 
   No modifications to this code may be distributed. If you do make any 
   modifications, please send them to SciTech Software so they may be 
   included in future releases.
   
   Permission to distribute this code as part of other programs may be 
   obtained on request from SciTech Software.
   
**************************************************************************

   Description:
   ============
   Hard formatted version of sscanf(). Implements FORTRAN-like file 
   reading.
 
   The only parsing characters recognised are:
      %<n>f    A single precision floating point number of width <n>
      %<n>lf   A double precision floating point number of width <n>
      %<n>d    An integer of width <n>
      %<n>ld   A long integer of width <n>
      %<n>u    An unsigned of width <n>
      %<n>lu   An unsigned long of width <n>
      %<n>s    A string of width <n>
      %c       A character (of width 1)
      %<n>x    <n> spaces (like FORTRAN).
   With the exception of the %c parser, the column width, <n>,
   *must* be specified.

   Blank fields read as numbers are given a value of zero.

   Returns: The number of arguments filled in (0 if end of file or no
            specifiers found in format string).



**************************************************************************

   Usage:
   ======
   For example:

   double   DoubVar;
   int      IntVar;
   char     CharVar,
            StringVar[16];

   fsscanf(buffer,"%8lf%5x%3d%c%3x%8s",
           &DoubVar,&IntVar,&CharVar,StringVar);


**************************************************************************

   Revision History:
   =================
   V1.0  17.06.93 Original    By: ACRM
   V1.1  12.07.93 Added %u and %lu. Corrected %s and %c to blank rather 
                  than NULL strings if buffer runs out. Pads string if 
                  buffer ran out in the middle. Takes \n in buffer as end 
                  of string.

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <ctype.h>

#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/*>int fsscanf(FILE *fp, char *format, ...)
   ----------------------------------------
   Hard formatted version of fscanf(). Implements FORTRAN-like file 
   reading.
 
   The only parsing characters recognised are:
      %<n>f    A single precision floating point number of width <n>
      %<n>lf   A double precision floating point number of width <n>
      %<n>d    An integer of width <n>
      %<n>ld   A long integer of width <n>
      %<n>u    An unsigned of width <n>
      %<n>lu   An unsigned long of width <n>
      %<n>s    A string of width <n>
      %c       A character (of width 1)
      %<n>x    <n> spaces (like FORTRAN).
   With the exception of the %c parser, the column width, <n>,
   *must* be specified.

   Blank fields read as numbers are given a value of zero.

   Returns: The number of arguments filled in (0 if end of file or no
            specifiers found in format string).

   17.06.93 Original    By: ACRM
   12.07.93 Added %u and %lu. Corrected %s and %c to blank rather than
            NULL strings if buffer runs out. Pads string if buffer ran
            out in the middle. Takes \n in buffer as end of string.
*/
int fsscanf(char *buffer, char *format, ...)
{
   va_list        ap;
   char           *FormStart,
                  *BuffStart,
                  *stop,
                  form[16],
                  value[40],
                  *ptr,
                  type;
   int            i,
                  *IntPtr,
                  NArg     = 0,
                  width    = 0;
   BOOL           LongType = FALSE;
   double         *DblPtr;
   float          *FloatPtr;
   long           *LongPtr;
   unsigned       *UPtr;
   unsigned long  *ULongPtr;

   /* Blank the buffers                                                 */
   memset(value,  '\0', 40);

   /* Return if line is blank                                           */
   if(!buffer || buffer[0] == '\0' || buffer[0] == '\n') return(EOF);
   
   /* Start the variable argument processing                            */
   va_start(ap, format);

   /* Intialise FormStart to the start of the format string and BuffStart
      to start of input buffer
   */
   FormStart = format;
   BuffStart = buffer;

   for(;;)
   {
      /* Flag for long variables                                        */
      LongType = FALSE;
      
      /* Find the start of a % group from the format string             */
      while(*FormStart && *FormStart != '%') FormStart++;
      if(!(*FormStart)) break;      /* Exit routine                     */
   
      /* Find the next occurence of a %                                 */
      stop = FormStart+1;
      while(*stop && *stop != '%') stop++;
   
      /* Copy these format characters into our working buffer           */
      for(i=0; FormStart != stop; i++)
         form[i] = *(FormStart++);
      form[i] = '\0';

      /* Find the type we're dealing with                               */
      ptr = form + i;
      while(*ptr == '\0' || *ptr == ' ' || *ptr == '\t') ptr--;
      type = toupper(*ptr);
      
      /* Set long flag if appropriate                                   */
      if((*(ptr-1) == 'l') || (*(ptr-1) == 'L'))
         LongType = TRUE;

      /* If it's not a character, read the width from the form string   */
      width = 0;
      if(type == 'C')
      {
         width = 1;
      }
      else
      {
         for(ptr = form+1; *ptr && isdigit(*ptr); ptr++)
         {
            width *= 10;
            width += (*ptr) - '0';
         }
      }
      
      /* Extract width characters from the input buffer. If the input 
         buffer has run out, value will be a NULL string.
      */
      stop = BuffStart + width;
      for(i=0; *BuffStart && *BuffStart != '\n' && BuffStart != stop; i++)
         value[i] = *(BuffStart++);
      value[i] = '\0';
      
      /* Act on each type                                               */
      switch(type)
      {
      case 'F':      /* A double precision or float                     */
         if(LongType)
         {
            DblPtr = va_arg(ap, double *);
            if(sscanf(value,"%lf", DblPtr) == (-1))
               *DblPtr = (double)0.0;
         }
         else
         {
            FloatPtr  = va_arg(ap, float  *);
            if(sscanf(value,"%f",  FloatPtr) == (-1))
               *FloatPtr = (float)0.0;
         }

         break;
      case 'D':      /* An integer or long int                          */
         if(LongType)
         {
            LongPtr = va_arg(ap, long *);
            if(sscanf(value,"%ld", LongPtr) == (-1))
               *LongPtr = 0L;
         }
         else
         {
            IntPtr  = va_arg(ap, int  *);
            if(sscanf(value,"%d",  IntPtr) == (-1))
               *IntPtr = 0;
         }
         break;
      case 'U':      /* An unsigned or unsigned long                    */
         if(LongType)
         {
            ULongPtr = va_arg(ap, unsigned long *);
            if(sscanf(value,"%lu", ULongPtr) == (-1))
               *ULongPtr = 0L;
         }
         else
         {
            UPtr  = va_arg(ap, unsigned  *);
            if(sscanf(value,"%u",  UPtr) == (-1))
               *UPtr = 0;
         }
         break;
      case 'S':      /* A string                                        */
         ptr = va_arg(ap, char *);
         if(value[0])                  /* Input buffer not empty        */
         {
            *(value + width) = '\0';
            strncpy(ptr, value, width+1);

            /* If the input buffer ran out in this string, pad with 
               spaces and terminate.
            */
            if(strlen(ptr) < width) padterm(ptr, width);
         }
         else                          /* Input buffer empty            */
         {
            for(i=0; i<width; i++)
               *(ptr + i)  = ' ';
            *(ptr + width) = '\0';
         }
         break;
      case 'C':      /* A character (insert a space if buffer empty)    */
         *(va_arg(ap, char *)) = (value[0] ? value[0]: ' ');
         break;
      case 'X':      /* A column to skip                                */
         /* Fall through to default action                              */
      default:
         /* Do nothing                                                  */
         ;
      }
      
      /* If not a blank column, increment arg count                     */
      if(type != 'X') NArg++;
   }
   
   /* End variable argument parsing                                     */
   va_end(ap);

   /* Return number of values read                                      */
   return(NArg);
}

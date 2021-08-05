#include "ProtoTypes.h"
#include "CongenProto.h"
 
/*************************************************************************
parse()
Input:      char  *comline       A command line string to parse
            int   nkeys          Number of keywords
            KeyWd *keywords      Array of keyword structures
Output:     float *floatparam    Array of returned strings
            char  **strparam     Array of pointers to returned strings
Returns:    int                  Index of found command or error flag
*************************************************************************/
 
int parse(
char *comline,
int nkeys,
KeyWd *keywords,
float *floatparam,
char *strparam[]
)
{
   char *command;
   int  i,n,found,nletters,nlett;
 
   command = killspcs(comline);
   found = 0;
   if((command[0]=='!') ||
      (command[0]==LF)  ||
      (command[0]==CR)  ||
      (command[0]=='\0'))
      return(PARSE_COMMENT);
 
   for(i=0;i<nkeys;i++)
   {
      /* match() returns 1 if first string finishes first or exact match
                         2 if second string finishes first
                         0 if a mismatch
         We only want to act in the first case
      */
      if((n=match(command,(keywords[i]).name,&nletters))==1)
      {
#ifdef DEBUG
         printf("Matched letters: %d\n",nletters);
#endif
         if(found)  /* If found already */
         {
#ifdef NOISY
            fprintf(OUTPUT,"Ambiguous keyword: %s\n",command);
#endif
            return(PARSE_ERRC);
         }
         found = i+1;  /* +1, so keyword 0 will flag TRUE */
         nlett = nletters;
      }
   }
   if(!found)
   {
#ifdef NOISY
      fprintf(OUTPUT,"Keyword not known: %s\n",command);
#endif
      return(PARSE_ERRC);
   }
   command+=nlett;
   found--; /* Reset to point to the correct keyword */
 
#ifdef DEBUG
   printf("After getting command, line contains %s\n",command);
#endif
 
   /* Get data requirements for this keyword */
   if((keywords[found]).string)
   {
#ifdef DEBUG
      printf("Command expects strings\n");
#endif
      for(i=0; i<(keywords[found]).nparam; i++)
      {
         command = killspcs(command);
         if((nletters = GetString(command,strparam[i]))==0)
         {
#ifdef NOISY
            fprintf(OUTPUT,"Missing string parameter in: %s\n",comline);
            fprintf(OUTPUT,"Command expects %d parameters\n",
                   (keywords[found]).nparam);
#endif
            return(PARSE_ERRP);
         }
         command += nletters;
      }  /* End of for(i) */
   }
   else
   {
      /* A numeric or no parameter */
#ifdef DEBUG
      printf("Command expects %d parameters\n",(keywords[found]).nparam);
#endif
      for(i=0; i<(keywords[found]).nparam; i++)
      {
         command = killspcs(command);
         if(!GetParam(command,&(floatparam[i]),&nletters))
         {
#ifdef NOISY
            fprintf(OUTPUT,"Error in parameter %s\n",comline);
            fprintf(OUTPUT,"Command expects %d numeric parameters\n",
                   (keywords[found]).nparam);
#endif
            return(PARSE_ERRP);
         }
         command += nletters;
      }  /* End of for(i) */
   }  /* End of else */
   return(found);
}
 

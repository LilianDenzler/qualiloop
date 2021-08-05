#include "ProtoTypes.h"
#include "CongenProto.h"
 
/* Routine to convert 1-letter amino acid code to 3-letter code.
   Will set the 3-letter code to UNK if an unknown 1-letter code is
   specified.
 
   Input:  char *one        The one letter code
   Output: char *three      The three letter code
*/
 
void onethr(
char one,
char three[]
)
{
char  tab3[20][4],tab1[20];
int   i, found = 0;
 
   strcpy(tab3[0],"ALA");   tab1[0] = 'A';
   strcpy(tab3[1],"CYS");   tab1[1] = 'C';
   strcpy(tab3[2],"ASP");   tab1[2] = 'D';
   strcpy(tab3[3],"GLU");   tab1[3] = 'E';
   strcpy(tab3[4],"PHE");   tab1[4] = 'F';
   strcpy(tab3[5],"GLY");   tab1[5] = 'G';
   strcpy(tab3[6],"HIS");   tab1[6] = 'H';
   strcpy(tab3[7],"ILE");   tab1[7] = 'I';
   strcpy(tab3[8],"LYS");   tab1[8] = 'K';
   strcpy(tab3[9],"LEU");   tab1[9] = 'L';
   strcpy(tab3[10],"MET");  tab1[10] = 'M';
   strcpy(tab3[11],"ASN");  tab1[11] = 'N';
   strcpy(tab3[12],"PRO");  tab1[12] = 'P';
   strcpy(tab3[13],"GLN");  tab1[13] = 'Q';
   strcpy(tab3[14],"ARG");  tab1[14] = 'R';
   strcpy(tab3[15],"SER");  tab1[15] = 'S';
   strcpy(tab3[16],"THR");  tab1[16] = 'T';
   strcpy(tab3[17],"VAL");  tab1[17] = 'V';
   strcpy(tab3[18],"TRP");  tab1[18] = 'W';
   strcpy(tab3[19],"TYR");  tab1[19] = 'Y';
 
   for(i=0; i<20; i++)
   {
      if(tab1[i]==one)
      {
         strcpy(three,tab3[i]);
         found = 1;
         break;
      }
   }
   if(!found)
      strcpy(three,"UNK");
}
 

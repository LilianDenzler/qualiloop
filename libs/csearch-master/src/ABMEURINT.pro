#ifndef __ABMEURINT_PRO__
#define __ABMEURINT_PRO__
 
/*
   Prototypes for directory : ABMEURINT
   Generated by pimmsccs on Tue Nov  3 13:44:21 GMT 1992
   OMUPD rkw 02/11/93 Added aemain prototype
*/
 
#ifdef __FORT_UNDERSCORE__
#define aemain aemain_
#endif

void DoABMEURINT (int argc ,char * argv [ ]);

void aemain(char *cinpdb,int *linpdb,char *cinign,int *linign,char *cintyp,int *lintyp,char *coudef,int *loudef,char *coucri,int *loucri,char *coucon,int *loucon,char *cuoerr,int *luoerr);
 
#endif

#ifndef __SEQU_PRO__
#define __SEQU_PRO__

/***********************************************************************
 *      Name: SeqU.pro                                                 *
 *  Function: Declares prototypes for sequence alignment program       *
 *            utilities functions (housekeeping routines)              *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: KY Cockwell, Oxford Molecular Ltd.                       *
 *      Date: 05/07/19931                                              *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 * 16/09/93 DW      Changed prototype on SeqUfreeALN()                 *
 * 23/09/93 KYC     Add prototypes for Secondary structure routines    *
 ***********************************************************************/

#include "SeqU.h"
#include "SeqCalc.h"

ALNPtr  SeqUallocALN( void);
ALNPtr  SeqUfreeALN( ALNPtr align);
SEQPtr  SeqUallocSEQ( void);
void    SeqUfreeSEQ( SEQPtr seq);
RESPtr  SeqUallocRES( void);
void    SeqUfreeRES( RESPtr res);
CGFPtr  SeqUallocCGF( void);
CGFPtr  SeqUfreeCGF( CGFPtr res);
CGFGPtr SeqUallocGrp( void);
void    SeqUfreeGrp( CGFGPtr res);
PRFPtr  SeqUallocPRF( void);
PRFPtr  SeqUfreePRF( PRFPtr prf);
MATPtr  SeqUallocMAT( void);
MATPtr  SeqUfreeMAT( MATPtr mat);
ProfPtr SeqUallocPfl( void);
ProfPtr SeqUfreePfl( ProfPtr profile);
MOTPtr  SeqUallocMOT( void);
MOTPtr  SeqUfreeMOT( MOTPtr mot);
SecSPtr SeqUallocSS(SEQPtr seq);
SecSPtr SeqUfreeSS(SEQPtr seq, SecSPtr oldptr);

SEQPtr  SeqUfindACC( ALNPtr align, char *accession);
SEQPtr  SeqUfindN( ALNPtr align, short n);
RESPtr  SeqUfindRES( SEQPtr seq, short n);
RESPtr  SeqUlastRES( SEQPtr seq);
void    SeqUdiscSEQ( ALNPtr align, SEQPtr seq);
SEQPtr  SeqUdiscN( ALNPtr align,short pos);

short   SeqUcompACC( char *access1,char *access2);

void    SeqUinsSEQ( ALNPtr align, SEQPtr seq, short pos);
void    SeqUappSEQ( ALNPtr align, SEQPtr seq);

void    SeqUinsRES( SEQPtr seq, RESPtr res, short pos);
void    SeqUappRES( SEQPtr seq, RESPtr res);
void    SeqUdelRES( RESPtr res, SEQPtr seq);
void    SeqUdelRESN( SEQPtr seq, short pos);
void    SeqUsetRES( RESPtr res, char amino);

short   SeqUres2idx( char amino);
char    SeqUidx2res( short index);

char   *SeqUsequence( SEQPtr seq, short gap);
char   *SeqUsegment( SEQPtr seq, short start, short length, short gap);
SEQPtr SeqUstrSEQ(char *text, char *code, char *organism);

void    SeqUgroup( SEQPtr seq, short group);
void    SeqUungroup( SEQPtr seq);
void    SeqUsetflag( SEQPtr seq, short flag, short state);
short   SeqUgetflag( SEQPtr seq, short flag);
void    SeqUselect( SEQPtr seq);
void    SeqUallsel( ALNPtr align);
void    SeqUallclear( ALNPtr align);
void    SeqUdeselect( SEQPtr seq);
void    SeqUlock( SEQPtr seq);
void    SeqUunlock( SEQPtr seq);
void    SeqUallunlk( ALNPtr align);

short   SeqUnSEQ( ALNPtr align);
short   SeqUnSel( ALNPtr align);
short   SeqUmaxSEQ( ALNPtr align, short selected);

char   *SeqUcon( ALNPtr align);
char   *SeqUconSEG( SEQPtr conseq, short start, short length);
char   *SeqUident( ALNPtr align);

void    SeqUcalccon( ALNPtr align);
void    SeqUcalcid( ALNPtr align);
void    SeqUconsel( ALNPtr align, SEQPtr consens, SEQPtr ident);

CGFGPtr SeqUfindCGF( CGFPtr cgf, char amino);

short   SeqUfillGrp( CGFGPtr group, char *comment, char *rgbamino);
short   SeqUscore( MATPtr mat, char amino1, char amino2);

float   SeqUfindPRF( PRFPtr prf, char amino);

ProfPtr SeqUcalcPRF( SEQPtr seq, PRFPtr prf, short L);
void    SeqUprofcalc( ProfPtr profile, PRFPtr prf, SEQPtr seq);
void    SeqUrclcProf( ProfPtr profile, short L);
ProfPtr SeqUfindProf( SEQPtr seq, char *name);
short   SeqUlookProf( ProfPtr profile, short n, float *value);

void    SeqUPrfStrp( SEQPtr seq);
void    SeqUgapSEQ( SEQPtr seq);
ALNPtr  SeqUresetALN(ALNPtr align);

ALNPtr  SeqUcopyALN(ALNPtr source);
SEQPtr  SeqUcopySEQ(SEQPtr source, short start, short end );
CGFPtr  SeqUcopyCGF(CGFPtr source);
LockPtr SeqUcopylck(LockPtr source, short nlocks);
CGFGPtr SeqUcopyGrp(CGFGPtr source);
short   SeqUappALN(ALNPtr source, ALNPtr target, short start, short end);
short   SeqUstemSEQ(ALNPtr source, ALNPtr target);

void    SeqUdelSSres( SEQPtr seq, short pos);
void    SeqUdelprofN( SEQPtr seq, short pos, short nres);

void    SeqUinsSSN(SEQPtr seq, short pos);
void    SeqUinsprofN(SEQPtr seq, short pos);
void    SeqUdeselall(ALNPtr align);
short   SeqUfindSPos(ALNPtr aln, SEQPtr seq);

void    SeqUsetptype(short cutoff, ALNPtr align);
void    SeqUsetpres(short cutoff, ALNPtr align);
void    SeqUdefRgn(SEQPtr seq, short rgn, short end, short r, short g, short b);
void    SeqUsetRgn(SEQPtr seq, int nreg, ...);
short   SeqUresNum(RESPtr rS, short gap);
char   *SeqUstrdup(char *original);
short   SeqUmaxRes(ALNPtr align);
short   SeqUmaxRng(ALNPtr align, short start, short end);
void    SeqUclnAcc(SEQPtr seq);

#endif

#ifndef __SEQIOSEQ__
#define __SEQIOSEQ__

/***********************************************************************
 *      Name: SeqIOseq.pro                                             *
 *  Function: Declares prototypes for sequence alignment program       *
 *            file io functions                                        *
 * Copyright: (C) Oxford Molecular Limited                             *
 *---------------------------------------------------------------------*
 *    Author: KY Cockwell, Oxford Molecular Ltd.                       *
 *      Date: 05/07/1993                                               *
 *---------------------------------------------------------------------*
 * Modification Record                                                 *
 * Date     Inits   Comments                                           *
 * 16/09/93 DW      Include SeqIOseq.h                                 *
 * 26/11/93 DW      Sequence/alignment I/O review (CR 2757)            *
 ***********************************************************************/

#include <stdio.h> 
#include <stdlib.h>
#include "SeqIOseq.h"
#include "SeqU.h"

SEQPtr SeqIOreadSEQ( char *filen, short *error);
ALNPtr SeqIOreadALN( char *filen, short type, short *error);
ALNPtr SeqIOrdMult( char *filen, short *error);
ALNPtr SeqIOrdSomap( char *filen, short *error);
ALNPtr SeqIOrdGCG( char *filen, short *error);
ALNPtr SeqIOrdClust(char *filen, short *error);
SEQPtr SeqIOrdSWISS( FILE *seqfile, short *error);
SEQPtr SeqIOreadPIR( FILE *seqfile, short *error);

short SeqIOwrtPIR(SEQPtr seq, char *filename);
short  SeqIOwrtALN( ALNPtr aln, char *file, short type,
                    short over, short pagesize, short solid);
short  SeqIOwrtMult( ALNPtr aln, char *file, short over);
short  SeqIOwrtASC( ALNPtr aln, char *file, short over, short pagesize, 
                    SEQPtr consens, SEQPtr ident);
short  SeqIOwrtRUL( FILE *outfil, short startpos, short pagesize,
                    short startcount);
short  SeqIOwrhead( FILE *outfil, short nseqs, ALNPtr align);
short  SeqIOwrtGCG( ALNPtr aln, char *file, short over);
short  SeqIOrdaln1(FILE *fPtr, char *buffer, ALNPtr newA, short sPos,
                   short aPos, short aEnd);
short  SeqIOrdaln2(FILE *fPtr, char *buffer, ALNPtr newA, short sPos,
                    short aPos, short aEnd);
char  *SeqIOgetacc(char *buffer, short aPos);
void   SeqIOtailSEQ( ALNPtr align);

FILE  *SeqIOopnfil( char *filename, char *mode, short over);
void   SeqIOclosef( FILE *filptr);

CGFPtr SeqIOrdCGF( char *filen, short *error);
char  *SeqIOgetR1( short pagesize, short startcount,
                   short startpos, short *error);
char  *SeqIOgetR2( short pagesize, short startpos, 
                   short startcount,  short *error);

short  SeqIOgetSEQ( SEQPtr currentseq, short pagesize,char *outline,
                    char *seqfin, RESPtr *reslist, short *finished,
                    short *seqcount);
short  SeqIOinitALN( ALNPtr align, char **seqfin, RESPtr **reslist, 
                     short *error);

short  SeqIOwrtPS(ALNPtr align, char *filename, short target, SEQPtr consens,
                  SEQPtr ident, short state);
short  SeqIOinitPS( char *filename, short target, int *Ysofar, int *Height);
short  SeqIOpage1( int *Ysofar, int Height, short *pagen, short nseqs);
short  SeqIOpage2( int *Ysofar, int Height, short *pagen, short nseqs);
short  SeqIOpsruler( short pagesize, int *Ysofar, short rulercount);
short  SeqIOpstail(ALNPtr align, short *idpos, int *Ysofar, short pagesize,
                   SEQPtr consens, SEQPtr ident, short state);
short  SeqIOpsRES( short seqcount, short pagesize, int Ysofar,
                   short *finished, char *seqfin, RESPtr *reslist, short nseqs,
                   short *EndAlign, CGFPtr cgf, short state );

SEQPtr SeqIOmultacc( FILE *alnfile, short nseqs, short *eCode);
short  SeqIOmltn( char *buffer);

MATPtr SeqIOrdMAT( char *filen, short *error);
void   SeqIOMATln( char *buffer, MATPtr newmat, short *aacodes, short ycount);
PRFPtr SeqIOrdPRF( char *file, short *error);

void   SeqIOselProf( SEQPtr seq, char *name, short state);
void   SeqIOselAllP( SEQPtr seq, short state);
short  SeqIOwrtProf( SEQPtr seq, char *title, char *file, short gap);

short  SeqIOcolCON(RESPtr res, CGFPtr cgf, short solid);
short  SeqIOcolPS(RESPtr res, CGFPtr cgf);
short  SeqIOfgets(char *buffer, int bufSize, FILE *fPtr);
RESPtr SeqIOprocRES(char resC, short *eCode);
SEQPtr SeqIOswissHdr(FILE *fPtr, char *buffer, short *eCode);
short  SeqIOrdSRes(FILE *fPtr, SEQPtr newS, char *endC);
short  SeqIOmltBlk(FILE *fPtr, char *buffer, ALNPtr newA);
short  SeqIOaccOK(char *buffer, short start, short end);

#endif

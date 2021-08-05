#ifndef __ABM_SEARCHDB_PRO__
#define __ABM_SEARCHDB_PRO__
/*****************************************************************************
 *      Name: SearchDB.pro                                                   *
 *  Function: Function prototypes for the AbM database search utility.       *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 20/04/93                                                       *
 *---------------------------------------------------------------------------*
 *    Inputs:                                                                *
 *   Outputs:                                                                *
 *   Returns:                                                                *
 * Externals:                                                                *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/
 
void SDBAddHit(CadEntry *TopEntry);
 
int SDBCheckLast(CadEntry *LastEntry);
 
int SDBCheckTop(CadEntry *TopEntry, int *ChainTooShort);
 
void SDBChkCAdis(CadHit **CurHit);
 
void SDBCopyResInf(int iamino);
 
void SDBCreateHit(CadHit **TopHitPtr, CadHit **LastHitPtr);
 
int SDBDoClust2(float Resolution);
 
void SDBGeneLoop(void);
 
void SDBGetAngChg(int FirstCA, int CompCA, float TarDis, int ActN, int ActCA, int ActC);
 
int SDBGetConst(char *CBuffer);
 
void SDBGetLoopCG(float LoopCG[3] );
 
int SDBGetPDBline(FILE *FilePtr, char *AtmNam, int *ResNum, float Coords[3]);
 
float SDBGetTor(float ix, float iy, float iz, float jx, float jy, float jz, float kx, float ky, float kz, float lx, float ly, float lz);
 
int SDBGetTorChg(int StatAt, int MoveAt, float TarDis, int AxAt1, int AxAt2, float RotAx[3], float *RotAng);
 
void SDBGetTrigVec(float X1, float X2, float X3, float Y1, float Y2, float Y3, float Z1, float Z2, float Z3, float YXAang, float XAlen, float XAvec[3]);
 
FILE *SDBInitCGfile(void);
 
void SDBMainProg(int argc, char *argv[]);
 
void SDBMakEntry(CadEntry **TopPtr, CadEntry **LastPtr);
 
void SDBMakeLoop(CadHit **CurHit);
 
FILE *SDBOpenFile(char *FileName, char *Access);
 
void SDBOrientRes(void);
 
int SDBQunique(CadHit *HitPtr1, CadHit *HitPtr2, float Toler);
 
int SDBRdCntrl(char *CntrlFile);
 
void SDBRdFrame(void);
 
void SDBRdResInfo(void);
 
int SDBReadEntry(FILE *CADFilePtr, CadEntry **TopPtr, CadEntry **LastPtr, long ExpRecLen);
 
long SDBRecLen(void);
 
void SDBRemAllEnt(CadEntry **TopPtr);
 
void SDBRemHit(CadHit **HitPtr);
 
void SDBRemHits(CadHit **TopHitPtr);
 
void SDBRemTopEnt(CadEntry **TopPtr);
 
void SDBRotFrag(float Xcart[], float Ycart[], float Zcart[], float Axis[3], float Angle, float RotCent[3], int Start, int Finish);
 
void SDBShowHits(void);
 
void SDBWrTors(void);
 
void SDBWriteLoop(FILE *CGfilePtr);
 
int SDBcluster(void);
 
void SDBcluster2(void);
 
void SDBfixCofG(void);
 
void SDBfixNCdis(void);
 
void SDBfixTerms(void);
 
void SDBoverlap(void);
 
void SDBsearch(void);
 
void SDBsuccess(void);

void SDBCrosProd( float Vec1[3], float Vec2[3], float Vec3[3]);

float SDBarccos(float CosAng);
 
#endif

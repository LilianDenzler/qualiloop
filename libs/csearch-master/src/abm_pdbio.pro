#ifndef __ABM_PDBIO_PRO__
#define __ABM_PDBIO_PRO__
/*****************************************************************************
 *      Name: abm_pdbio.pro                                                  *
 *  Function: Prototype definitions for the AbM PDB read/write routines      *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 15/10/93                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/
 
AtmData *PDBAllocAtm(AtmData **TopAtm, AtmData *AtmPtr, int AddCode);
 
ChnData *PDBAllocChn(ChnData **TopChn, ChnData *ChnPtr, int AddCode);
 
ResData *PDBAllocRes(ResData **TopRes, ResData *ResPtr, int AddCode);
 
int PDBChkAliSeq(PDBdata *PDBptr);
 
void PDBFreePDB(PDBdata *PDBptr);
 
char PDBGetAAcode(char *ResNam);
 
int PDBGetCoords(ResData *ResPtr, char *AtmNam, float *X, float *Y, float *Z);
 
int PDBProcAtom(PDBdata *PDBptr, char *CBuffer, int LineNo, char *ErrMsg);
 
int PDBProcLine(char *CBuffer, int *AtmNum, char *AtmNam, char *ResNam, int *ResNum, char *ChainID, char *Insert, float *Xcart, float *Ycart, float *Zcart, int LineNo, char *ErrMsg);
 
int PDBProcRem(PDBdata *PDBptr, char *CBuffer, int LineNo, char *ErrMsg);
 
PDBdata *PDBReadFile(char *FileName, char *ErrMsg);
 
void PDBRemAtm(PDBdata *PDBptr, char *AtmNam);
 
void PDBRemHyd(PDBdata *PDBptr);
 
void PDBRemRes(PDBdata *PDBptr, char *ResName);
 
void PDBRenAllAtm(PDBdata *PDBptr, char *ResNam, char *OldNam, char *NewNam);
 
void PDBRenResAtm(ResData *ResPtr, char *OldNam, char *NewNam);
 
void PDBRenumAtm(PDBdata *PDBptr);
 
void PDBRenumRes(PDBdata *PDBptr);
 
int PDBWriteFile(char *FileName, PDBdata *PDBptr, char *ErrMsg);
 
 
#endif

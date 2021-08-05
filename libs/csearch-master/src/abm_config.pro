#ifndef __ABM_CONFIG_PRO__
#define __ABM_CONFIG_PRO__
 
/*****************************************************************************
 *      Name: abm_config.pro                                                 *
 *  Function: Configuration file reading routine prototypes                  *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 14/06/93                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/

 
int DoReadConf(FILE *FilePtr, CONFIG *ConfPtr, char *ErrMess);
 
void GenExtraInf(CONFIG *ConfPtr);
 
CHAINCONF *GetChainPtr(CONFIG *ConfPtr, char ChainID);
 
FRAGCONF *GetFragPtr(CONFIG *ConfPtr, char *FragID);
 
int RdDataFiles(FILE *FilePtr, CONFIG *ConfPtr, int *LineNo, char *ErrMess);
 
int ReadChain(FILE *FilePtr, CHAINCONF *ChainPtr, int *LineNo, char *ErrMess);
 
CONFIG *ReadConfig(char *FileName, char *ErrMess);
 
int ReadFrags(FILE *FilePtr, FRAGCONF *FragPtr, int *LineNo, char *ErrMess);
 
#endif

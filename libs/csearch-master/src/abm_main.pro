#ifndef __abm_main_PRO__
#define __abm_main_PRO__
 
/*****************************************************************************
 *      Name: abm_main.pro                                                   *
 *  Function: The prototypes for functions in ABM_MAIN and ABM_MAIN2 packages*
 * Copyright: (C) OXFORD MOLECULAR LTD, 1993                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 23/11/93                                                       *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/

int ABMActOnCom(void);
 
int ABMAliToSeq(PDBdata *PDBptr);
 
void ABMAllUpper(char *TextString);
 
MenuData *ABMAllocMenu(int Ngroups, int Nareas[]);
 
IntrData *ABMAllocRes(LoopData *LoopPtr);
 
int ABMAskBldMeth(LoopData *LoopPtr);
 
int ABMAskBldOrd(LoopData *CurLoop);
 
void ABMBuildItems(short NumCmd, char *ArgValue);
 
int ABMBuildMenu(char *LineInput);
 
void ABMCGENItems(short NumCmd, char *ArgValue);
 
int ABMCGENMenu(char *LineInput);
 
void ABMChkAdvice(char *Name, int Start, int Finish, int Length, char *Method, char *Advice);
 
int ABMConsMenu(void);
 
int ABMCountRes(int FstRes, int LstRes, int ChainNum);
 
int ABMDBSrchMenu(char *LineInput);
 
void ABMDbaseItems(short NumCmd, char *ArgValue);
 
int ABMDbaseMenu(char *LineInput);
 
void ABMDoAdvice(void);
 
void ABMDrawMenu(MenuData *MenuPtr);
 
void ABMDrawRes(IntrData *ResPtr, LoopData *LoopPtr);
 
int ABMEnergyMenu(char *LineInput);
 
void ABMEngyItems(short NumCmd, char *ArgValue);
 
int ABMExtraInput(char *ArgPrompt, char *ArgValue);
 
int ABMFileExists(char *FileName);
 
void ABMFileItems(short NumCmd, char *ArgValue);
 
int ABMFileMenu(char *LineInput);
 
void ABMFinishUDF(void);
 
void ABMFireShell(void);
 
int ABMFrParmMenu(char *LineInput);
 
void ABMFrameItems(short NumCmd, char *ArgValue);
 
void ABMFreeMenu(MenuData *MenuPtr);
 
void ABMGetCGENres(LoopData *LoopPtr, int *Rebuild1, int *Rebuild2, int *Closure1, int *Closure2);
 
void ABMGetCanons(void);
 
int ABMGetFileName(char *ArgValue, int Qexist, int QOwrite,  int QAppend, char *FileType);
 
int ABMGetSystem(int QSingleWin);
 
void ABMGetTitle(char *LeftText, char *RiteText, char *TitleText);
 
int ABMGetUDBCalf(void);
 
int ABMGetUDFLine(FILE *FilePtr, char *KeyWord, char *ArgValue, int *LineNo);
 
void ABMHiLiteLoop(MenuData *MenuPtr);
 
void ABMHiLiteRes(IntrData *ResPtr, int QHiLite);
 
int ABMInclMenu(void);
 
int ABMInitialise(void);
 
int ABMIsRunning(void);
 
void ABMKeyWords(void);
 
void ABMLoopItems(short NumCmd, char *ArgValue);
 
void ABMLpDBItems(short NumCmd, char *ArgValue);
 
int ABMLpParmMenu(char *LineInput);
 
void ABMMPskel(void);
 
void ABMMainItems(short NumCmd, char *ArgValue);
 
int ABMMainMenu(char *LineInput);
 
void ABMMainProg(int argc, char *argv[]);
 
void ABMMakeDB(void);
 
int ABMMtchFrame(int ChainNum, char *UDBname);
 
LoopData *ABMMtchLoop(char *LoopNam);
 
void ABMNullSeqs(void);
 
int ABMPDBtoSeq(PDBdata *PDBptr);
 
char *ABMPirToAli(char *FileName);
 
void ABMProcess(void);
 
int ABMQaligned(void);
 
int ABMRdChnUDF(FILE *FilePtr, int *LineNo);
 
int ABMRdConsUDF(char *KeyWord, char *ArgList, LoopData *LoopPtr);
 
int ABMRdGenUDF(FILE *FilePtr, int *LineNo);
 
int ABMRdIncUDF(FILE *FilePtr, int ChnNum, int *LineNo, LoopData *LoopPtr);
 
int ABMRdLoopUDF(LoopData *LoopPtr, FILE *FilePtr, int *LineNo);
 
int ABMRdSChnUDF(FILE *FilePtr, int ChnNum, int *LineNo, LoopData *LoopPtr);
 
void ABMReModel(char *FileName);
 
int ABMReadConfig(void);
 
void ABMReadPDB(char *FileName);
 
void ABMReadSeq(char *FileName);
 
int ABMReadUDB(void);
 
void ABMReadUDF(char *FileName);
 
void ABMResCurLR(IntrData *ResPtr, LoopData *LoopPtr, int LorR, int MultiStep);
 
void ABMResCurUD(IntrData *ResPtr, LoopData *LoopPtr, int UorD);
 
int ABMResScrn(IntrData *ResPtr, LoopData *LoopPtr);
 
void ABMResetSys(void);
 
void ABMRunBuilder(int QReStart, int QWrControl);
 
int ABMSChnMenu(void);
 
void ABMSaveUDF(char *FileName);
 
void ABMSetAllCan(LoopData *LoopPtr);
 
MenuData *ABMSetBldMenu(void);
 
void ABMSetBldSeq(LoopData *LoopPtr, char *SeqText);
 
void ABMSetBuild(LoopData *LoopPtr);
 
MenuData *ABMSetCGENmenu(LoopData *LoopPtr);
 
void ABMSetClosure(LoopData *LoopPtr);
 
MenuData *ABMSetConsMenu(LoopData *LoopPtr, char *MsgText);
 
MenuData *ABMSetDbaseMenu(void);
 
void ABMSetDefBld(int QSetFrame);
 
MenuData *ABMSetEngyMenu(LoopData *LoopPtr);
 
MenuData *ABMSetFileMenu(void);
 
MenuData *ABMSetFrParm(void);
 
int ABMSetHiLite(MenuData *MenuPtr);
 
void ABMSetInclude(LoopData *LoopPtr);
 
int ABMSetLoopDef(LoopData *LoopPtr);
 
MenuData *ABMSetLpDBmenu(LoopData *LoopPtr);
 
MenuData *ABMSetLpMenu(LoopData *LoopPtr);

void ABMSetLpMsg(MenuData *MenuPtr);
 
MenuData *ABMSetMainMenu(void);
 
void ABMSetOrder(void);
 
void ABMSetRBrange(LoopData *LoopPtr);
 
void ABMSetResMsg(IntrData *ResPtr);
 
void ABMSetResVal(int ValueFlag, IntrData *ResPtr, LoopData *LoopPtr);
 
void ABMSetSChains(LoopData *LoopPtr);
 
MenuData *ABMSetSysMenu(void);
 
void ABMSetUpRes(IntrData *ResPtr, LoopData *LoopPtr);
 
void ABMShStatus(void);
 
void ABMShutDown(void);
 
void ABMStepClosure(MenuData *MenuPtr, int StepMode);
 
void ABMStepRBrange(MenuData *MenuPtr, int StepMode);
 
void ABMSysItems(short NumCmd, char *ArgValue);
 
int ABMSysMenu(char *LineInput);
 
int ABMUserInput(MenuData *MenuPtr, int StrMode, char *Reply, int *ActGrp, int *ActArea);
 
int ABMValidRes(IntrData *ResPtr, LoopData *LoopPtr);
 
int ABMWrControl(char *FCntrl);
 
void ABMWrIgnRes(FILE *FilePtr, LoopData *LoopPtr);
 
void ABMWrLoopBld(FILE *FPtr, LoopData *LpPtr);
 
void ABMWrSideChn(FILE *FPtr, LoopData *LpPtr);
 
void ABMWrSideRes(FILE *FilePtr, LoopData *LoopPtr);
 
void ABMWriteMsg(char *ModName, char *MsgID, char *Embed1, char *Embed2, char *Embed3, char *Embed4);
 
void ABMWriteSeq(char *FileName);
 
void ABMcentre(char *InText, int TextLen, char *OutText);
 
int ABMconfirm(char *Message, int Qdefault);
 
 
#endif

int main(int argc, char **argv)
;
void SetDefaults(EPARAMS *EParams, FLAGS *Flags, char *ParamFile, 
                 char *TopFile)
;
BOOL OpenFiles(char *PDBFile, char *ControlFile, FILE **pdbfp, 
               FILE **controlfp)
;
BOOL ParseCmdLine(int argc, char **argv, char *PDBFile, char *ControlFile)
;
void UsageExit(void)
;
BOOL ParseControlFile(FILE *controlfp, EPARAMS *EParams, char *ParamFile,
                      char *TopFile, FLAGS *Flags, FILE **conffp, 
                      FILE **pdbfp, FILE **outfp)
;
ZONE *StoreZone(ZONE **zone, char *resspec1, char *resspec2)
;
BOOL BuildCache(int size)
;
void CacheConf(int ConfNum, REAL energy, int CacheSize)
;
void ShowCache(FILE *out, int NCache)
;
BOOL StoreUserDisulphide(char *resspec1, char *resspec2)
;

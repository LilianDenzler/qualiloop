TOPOLOGY *BuildTop(char *seq)
;
BOOL BuildAtoms(TOPOLOGY *topol, char *sequence)
;
int NAtomsInSequence(char *sequence)
;
BOOL BuildRest(TOPOLOGY *topol, char *sequence)
;
ATOM *FindAtomInRange(char *atom, TOPOLOGY *topol, int start, int stop)
;
char *OneThrTer(char one)
;
BOOL PatchTermini(TOPOLOGY *topol)
;
BOOL PatchTerminalNames(TOPOLOGY *topol)
;
BOOL PatchNTerDonor(TOPOLOGY *topol)
;
void DeleteTopology(TOPOLOGY *topol)
;
BOOL BuildBonds(TOPOLOGY *topol, char *resnam, int res, 
                int start, int stop)
;
BOOL BuildAngles(TOPOLOGY *topol, char *resnam, int res, 
                 int start, int stop)
;
BOOL BuildTorsions(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
;
BOOL BuildImpropers(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
;
BOOL BuildDonors(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
;
BOOL BuildAcceptors(TOPOLOGY *topol, char *resnam, int res, 
                   int start, int stop)
;
BOOL BuildExclusions(TOPOLOGY *topol, char *resnam, int res, 
                     int start, int stop)
;
TOPOLOGY *BuildDisulphide(TOPOLOGY *topol1, int res1, 
                          TOPOLOGY *topol2, int res2, 
                          MOLECULE *mol, TOPOLOGY *sstopol)
;
int ScanDisulphides(MOLECULE *mol)
;
BOOL PatchSSExclusions(TOPOLOGY *topol1, ATOM *AtomS1, ATOM *AtomCB1,
                       TOPOLOGY *topol2, ATOM *AtomS2, ATOM *AtomCB2)
;
int GetResFromAtomNum(TOPOLOGY *topol, int atomnum)
;
BOOL UserDisulphide(MOLECULE *mol, char *resspec1, char *resspec2)
;
BOOL PatchZones(MOLECULE *mol, ZONE *zones)
;
BOOL TrimTopology(MOLECULE *mol, ZONE *zones)
;
BOOL InZones(ATOM *atom, ZONE *zones)
;

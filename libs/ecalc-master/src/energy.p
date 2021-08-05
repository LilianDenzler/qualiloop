REAL EBond(TOPOLOGY *topol)
;
REAL EAngle(TOPOLOGY *topol)
;
REAL ETorsion(TOPOLOGY *topol)
;
REAL EImproper(TOPOLOGY *topol)
;
GRID *BuildGrid(MOLECULE *mol, REAL GridCut)
;
BOOL AddToGrid(GRID *grid, int AtomNum, ATOM *atom, REAL GridCutSq)
;
void TrimGrid(int NGrid, GRID *grid, ZONE *zones, ZONE *ignores)
;
REAL ENonBond(MOLECULE *mol, EPARAMS *eparams, 
              REAL *VdwAE, REAL *VdwRE, REAL *ElectE)
;
REAL EHBond(MOLECULE *mol, EPARAMS *eparams)
;
GRID *BuildHBGrid(MOLECULE *mol, REAL GridCut)
;
BOOL DoEnergyCalculations(FILE *out, MOLECULE *mol, EPARAMS EParams, 
                          FLAGS Flags, FILE *conffp)
;
REAL ShowEnergy(FILE *out, MOLECULE *mol, EPARAMS *eparams, FLAGS *flags)
;
BOOL GetConf(FILE *conffp, MOLECULE *mol, BOOL SkipReference)
;
ATOM **ReadConfHeader(FILE *conffp, MOLECULE *mol, BOOL SkipReference,
                      int *ncons)
;
REAL EResidue(ZONE *zone)
;
REAL ScoreResidue(char restype, REAL angle)
;
ATOM *FindDistalAtom(TOPOLOGY *topol, char restype, int start, int stop)
;
void FindSidechainCG(TOPOLOGY *topol, int start, int stop, VEC3F *CofG)
;
int AtomCount(ATOM **AtomArray, ATOM *atom)
;

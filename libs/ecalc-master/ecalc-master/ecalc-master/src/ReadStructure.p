MOLECULE *ReadStructure(FILE *fp, BOOL DoSS, BOOL *MissingAtoms, 
                        int *NumSS)
;
PDB *StoreCoords(TOPOLOGY *topol, PDB *pdb, BOOL *error)
;
BOOL CheckCoords(TOPOLOGY *topol)
;

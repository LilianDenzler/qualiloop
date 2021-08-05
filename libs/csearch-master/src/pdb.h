/* Include file for PDB handling
   *****************************
   Ammended 22.03.90: Added Secondary structure routines
   Ammended 22.03.90: Corrected field widths for V1.2 of ReadPDB()
   Ammended 22.03.90: Added clear_pdb()
   Ammended 22.03.90: Changed SEC structure to correct chain and ins widths
   Ammended 22.03.90: Added init_index macro
*/
struct pdb_entry
{
   char   junk[7];
   int    atnum;
   char   atnam[5];
   char   resnam[5];
   int    resnum;
   char   insert[2];
   char   chain[2];
   float  x,y,z,occ,bval;
   struct pdb_entry *next;
} ;
typedef struct pdb_entry PDB;

#define init_pdb(x) x = (PDB *)malloc(sizeof(PDB))

#define SELECT(x,w) x = (char *)malloc(5 * sizeof(char)); strcpy(x,w)

struct sec_entry
{
   char   chain1[2];
   int    res1;
   char   ins1[2];
   char   chain2[2];
   int    res2;
   char   ins2[2];
   char   type;
   struct sec_entry *next;
} ;
typedef struct sec_entry SEC;

#define init_sec(x) x = (SEC *)malloc(sizeof(SEC))

#define clear_pdb(p) strcpy(p->junk,"      "); \
                     p->atnum = 0; \
                     strcpy(p->atnam,"    "); \
                     strcpy(p->resnam,"    "); \
                     p->resnum = 0; \
                     strcpy(p->insert," "); \
                     strcpy(p->chain," "); \
                     p->x = 0.0; p->y = 0.0; p->z = 0.0; \
                     p->occ = 0.0; p->bval = 0.0; \
                     p->next = NULL

#define init_index(z,y) z=(PDB **)malloc(((y)+1) * sizeof(PDB *))

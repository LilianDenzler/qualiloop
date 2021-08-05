#include "ProtoTypes.h"
#include "CongenProto.h"
#include "values.h"
 
/***********************************************************************/
/* Read a PDB coordinate file into the coordinate arrays.
   Expects the coordinates to have NTER and CTER residues WITH atoms
   correctly specified.
   Handles ordering of atoms within residues.
   We use -9999.0 to specify atoms which have not had their coordinates
   set by the PDB file. This routine will exit if any atoms are left
   with coordinates like this. +9999.0 is used for undefined coordinates
   within CONGEN i.e. cleared coordinates which will be rebuilt.
*/
/***********************************************************************/
/* this flag is used to allow the use (or not) of adding H atoms etc to file */
extern int fixpdbfile;
 
void ReadCoords(
FILE *fp,
FILE *pgp_fp,
int showcoords
)
{
   PDB *pdb,
       *p;
   FILE *fpc;
   int natom,
       i, j;
 
   float x[20], y[20], z[20];
   char  resnam[8],
         chain[8],
         temp[8],
         atnam[20][4];
   int   resnum,
         rescount = 0,
         atcount,
         istart,
         istop;
 
 
   /* Read in the PDB file into a linked list */
   if(!fp)     prdie("ReadCoords() failed: PDB file not opened.\n");
   init_pdb(pdb);
   ReadPDB(fp,pdb,&natom);
 
   /* Add hydrogens and fix other coords */
   if(!pgp_fp) prdie("ReadCoords() failed: PGP file not opened.\n");
   if(fixpdbfile) fixpdb(pgp_fp,pdb);
 
   /* We don't need the PDB or PGP files any more, so close them */
   if(fp)
   {
      fclose(fp);
      fp = NULL;
   }
   if(pgp_fp)
   {
      fclose(pgp_fp);
      pgp_fp = NULL;
   }
 
 
   /* Set all coords to -9999.0 */
   for(i=0;i<maxat;i++)
      coords.xcart[i] = coords.ycart[i] = coords.zcart[i] = -9999.0;
 
 
   /* Step through the linked list */
   for(p=pdb;p;)
   {
      atcount=0;
 
      /* Copy the current residue */
      for(strcpy(chain,p->chain), strcpy(resnam,p->resnam), resnum = p->resnum;
          p && p->resnum == resnum && !strcmp(p->chain,chain);
          NEXT(p))
      {
         x[atcount] = p->x;
         y[atcount] = p->y;
         z[atcount] = p->z;
         strncpy(atnam[atcount],p->atnam,4);
         abmpad(atnam[atcount],4);
         atcount++;
      }
      rescount++;
 
      /* See if this is the type it's supposed to be according to the pstruct */
      padterm(resnam,4);
      if(strncmp(resnam,pstruct.resnme[rescount-1],4))
      {
         strncpy(temp,pstruct.resnme[rescount-1],4); temp[4] = '\0';
         fprintf(out,"Residue %d. Coordinate type %s, does not match topology \
type %s\n",rescount,resnam,temp);
         die();
      }
 
      /* Set the resid */
      sprintf(pstruct.resid[rescount-1],"%d",resnum);
      abmpad(pstruct.resid[rescount-1],4);
 
      /* Find the range of atoms spanned by this residue */
      istart = pstruct.lstatm[rescount-1];
      istop  = pstruct.lstatm[rescount];
 
      /* For each atom, get its coords by matching the atom name */
      for(i=istart; i<istop; i++)
      {
         for(j=0;j<atcount;j++)
         {
            if(!strncmp(atnam[j],pstruct.atmnme[i],4))
            {
               coords.xcart[i] = x[j];
               coords.ycart[i] = y[j];
               coords.zcart[i] = z[j];
               break;
            }
         }
      }
   }
 
   /* See if we've got any atoms left without defined coordinates */
   for(i=0;i<values.natoms;i++)
   {
      if(coords.xcart[i] == -9999.0)
      {
         fprintf(out,"Coordinates undefined for atom %d\n",i+1);
         for(p=pdb, j=0; p!=NULL && j<i; NEXT(p), j++);
         fprintf(out,"Error in residue %s %d, atom %s\n", 
                 p->resnam, p->resnum, p->atnam);
         die();
      }
   }
   
   if((fpc = fopen("reference.pdb","w+")) != NULL)
   {
      WritePDB(fpc,pdb);
      fclose(fpc);
   }   

   /* Free the pdb linked list as we've now finished with it */
   FREELIST(pdb,PDB);
}
 
 

#include <stdio.h>
#include "header.h"

struct {
       short int *space_grid;
       int ngridx,ngridy,ngridz;
       float xmn,ymn,zmn,xmx,ymx,zmx;
       float spgridsz,recipgrid;
       logical *excluded;
       int *cntnbx;
       logical *ingrid;
       int maxnbx;
       short int *nbxa;
       int freehd,freecls;
       int *nexthd,*clshd,*clstl,*nextcls,*clsatm;
       int *donp,*accp;
       int *resbya;
       logical *qside;
       float *radius,maxradius;
       } grid;

int main(int argc, char **argv)
{
   int i,j,k;
   
   grid.ngridx = 10;
   grid.ngridy = 10;
   grid.ngridz = 10;
   grid.space_grid = malloc(grid.ngridx*grid.ngridy*grid.ngridz*sizeof(short int));
   
   printf("space_grid malloc'd at: %d\n", grid.space_grid);
   fflush(stdout);

   for(i=0; i<grid.ngridx; i++)
   {
      for(j=0; j<grid.ngridy; j++)
      {
         for(k=0; k<grid.ngridz; k++)
         {
            int offset = (i * grid.ngridy * grid.ngridz) + (j * grid.ngridz) + k;
            
            grid.space_grid[offset] = 88;
         }
      }
   }
   
   test_for();

   printf("space_grid malloc'd at: %d\n", grid.space_grid);
   fflush(stdout);

   for(i=0; i<grid.ngridx; i++)
   {
      for(j=0; j<grid.ngridy; j++)
      {
         for(k=0; k<grid.ngridz; k++)
         {
            int offset = (i * grid.ngridy * grid.ngridz) + (j * grid.ngridz) + k;
            
            printf("grid[%d][%d][%d] = %d\n", i,j,k,grid.space_grid[offset]);
         }
      }
   }
   
      

   return(0);
   
}

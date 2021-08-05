#include "header.h"

struct 
{
   int n;
   short int *array;
}  grid;

int main(int argc, char **argv)
{
   int i = 99;
   
   grid.n = 20;
   grid.array = (short int *)malloc(grid.n * sizeof(short int));
   test_for(&i); 

   for(i=0; i<grid.n; i++)
      printf("grid.array[%d] = %d\n",i, grid.array[i]);
   
   return(0);
}

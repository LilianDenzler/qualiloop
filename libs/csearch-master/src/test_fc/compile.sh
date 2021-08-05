cc -g -c *.c
g77 -g -c *.f
cc -g -o mytest *.o -lm -lefence



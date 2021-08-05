#ifndef __RESTOPIO_H__
#define __RESTOPIO_H__


#ifdef MAIN

/* Global variables for ResTop reader */

char *ResTop_strparam[RESTOPMAXSTRPARAM]; /* Array for returned strings */
float ResTop_numparam[RESTOPMAXNUMPARAM];  /* Array for returned numbers */
int ResTop_errcnt;
int ResTop_wrncnt;
int nitc;

#else

/* Global variables for ResTop reader */

extern char *ResTop_strparam[RESTOPMAXSTRPARAM]; /* Array for returned strings */
extern float ResTop_numparam[RESTOPMAXNUMPARAM];  /* Array for returned numbers */
extern int ResTop_errcnt;
extern int ResTop_wrncnt;
extern int natc;
 
#endif

#endif

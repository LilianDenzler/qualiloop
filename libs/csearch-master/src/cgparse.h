#ifndef __CGPARSE_H__
#define __CGPARSE_H__

#define NSIDEOPT 5
#define SIDEOPTWIDTH 3
#define NEVALOPT 2
#define EVALOPTWIDTH 2

#ifdef MAIN
int nsideopt = NSIDEOPT;
int sideoptwidth = SIDEOPTWIDTH;
char st_sideopt[NSIDEOPT][4*SIDEOPTWIDTH] = {
           "FIRST      ",
           "INDEPENDENT",
           "ALL        ",
           "COMBINATION",
           "ITERATIVE  " };

int nevalopt = NEVALOPT;
int evaloptwidth = EVALOPTWIDTH;
char st_evalopt[NEVALOPT][4*EVALOPTWIDTH] = {
           "ENERGY ",
           "RMS    " };
int l_evalopt[NEVALOPT] = { 6,3 };

char *move_name[] = {
           "",
           "Chain Closure",
           "Backbone",
           "Side Chain",
           "Write Coordinates",
           "Status",
           "Read Best Conformations",
           "Evaluate" };

#else

extern int nsideopt;
extern int sideoptwidth;
extern char st_sideopt[NSIDEOPT][4*SIDEOPTWIDTH];
extern int nevalopt;
extern int evaloptwidth;
extern char st_evalopt[NEVALOPT][4*EVALOPTWIDTH];
extern int l_evalopt[NEVALOPT];
extern char *move_name[];
#endif

#endif

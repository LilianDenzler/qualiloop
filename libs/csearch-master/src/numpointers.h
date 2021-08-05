#ifndef __NUMPOINTERS_H__
#define __NUMPOINTERS_H__

#ifdef MAIN
int one   = 1,
    two   = 2,
    three = 3,
    four  = 4,
    five  = 5,
    six   = 6,
    seven = 7,
    eight = 8,
    nine  = 9,
    zero  = 0,
    eighty = 80;
float fzero = 0.0;
#else
extern int one,
           two,
           three,
           four,
           five,
           six,
           seven,
           eight,
           nine,
           zero,
           eighty;
extern float fzero;
#endif

#endif

#ifndef __CONGLOB_H__
#define __CONGLOB_H__

 
extern char   gres[MAXTYPE][8],
              gatom[MAXTYPE][8][8],
              grname[8],
              ghname[8][8],
              gnat[MAXATINRES][8];
extern int    igtype[MAXTYPE],
              info_level,
              npgp,
              no,
              kmax,
              it1, it2, it3, it4, it5, ih;
extern float  gr[MAXTYPE],
              alpha[MAXTYPE],
              beta[MAXTYPE],
              gx[MAXATINRES], gy[MAXATINRES], gz[MAXATINRES],
              hfac, fac;

#endif

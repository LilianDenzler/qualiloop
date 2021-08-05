#ifndef __GRID_H__
#define __GRID_H__

#define contact 1
#define nbond_energy 2
#define donor_hbond_energy 3
#define acceptor_hbond_energy 4
 
extern struct {
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

#endif

#include <stdio.h>
#include <stdlib.h>

#include "dtype.h"

/* graph name : temp_trans
//     4    5    6    7
//     o    |    |    |
//          |----|----|
//     o    |    |    |
//     0    1    2    3
//
// weight: 1
*/

int link_rule_temp_trans[]
        = {0,1,1,1,4,1,1,1,1,6,1,1,1,1,1,1};

int insert_rule_temp_trans(int* state) {
    if((state[1]==state[2])&&(state[2]==state[3])) {
        return 1;
    }
    return 0;
}

/* graph name : reference_conf
//     3    4    5 
//     |    |    |
//     |----|----|
//     |    |    |
//     0    1    2
//
// weight: -lambda
*/

int link_rule_reference_conf[]
        = {0,0,0,0,0,0,6,1,1,1,1,1};

int insert_rule_reference_conf(int* state) {
    if((state[0]==state[1])&&(state[1]==state[2])) {
        return 1;
    }
    return 0;
}

/* graph name : gauss_law
//     2    3
//     |    |
//     |----|
//     |    |
//     0    1
//
// weight: 0
*/

int link_rule_gauss_law[]
        = {0,0,0,0,4,1,1,1};

int insert_rule_gauss_law(int* state) {
    return 1;
}


static void create_cmf(double* cmf, double* weight, int length) {
    int i=0;

    cmf[0] = weight[i];
    for(i=1;i<length;i++) {
        cmf[i] = cmf[i-1]+weight[i];
    }
}

model* quantum_link_triangular_lattice(int lx, int ly, double lambda) {
    int nsite = 2*lx*ly;
    int nbond = 10*lx*ly;
    int mhnspin = 4;
    model* m = malloc_model(nsite,nbond,mhnspin);

    int n=0;
    // graph name : temp_trans
    // range [0*lx*ly, 1*lx*ly)
    // location   : sublattice B
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx          + x;                 //(x  ,y  ) sublattice A
            int j = y*lx          + (x+1)%lx;          //(x+1,y  ) sublattice A
            int k = ((y+1)%ly)*lx + (x+1)%lx;          //(x+1,y+1) sublattice A
            int l = y*lx          + (x+1)%lx + lx*ly;  //(x+1,y  ) sublattice B

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = 1.0;
            m->sweight += 1.0;
            m->bond2index[n*mhnspin+0] = l;
            m->bond2index[n*mhnspin+1] = i;
            m->bond2index[n*mhnspin+2] = j;
            m->bond2index[n*mhnspin+3] = k;
            n++;
        }
    }
    // graph name : temp_trans
    // range [1*lx*ly, 2*lx*ly)
    // location   : sublattice A
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx          + x + lx*ly;         //(x  ,y  ) sublattice B
            int j = ((y+1)%ly)*lx + x + lx*ly;         //(x  ,y+1) sublattice B
            int k = ((y+1)%ly)*lx + (x+1)%lx + lx*ly;  //(x+1,y+1) sublattice B
            int l = ((y+1)%ly)*lx + x;                 //(x  ,y+1) sublattice A

            m->bond2type[n]   = 0;
            m->bond2hNspin[n] = 4;
            m->bond2weight[n] = 1.0;
            m->sweight += 1.0;
            m->bond2index[n*mhnspin+0] = l;
            m->bond2index[n*mhnspin+1] = i;
            m->bond2index[n*mhnspin+2] = j;
            m->bond2index[n*mhnspin+3] = k;
            n++;
        }
    }
    // graph name : reference_vonf
    // range [2*lx*ly, 3*lx*ly)
    // location   : sublattice A
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx          + x;                 //(x  ,y  ) sublattice A
            int j = y*lx          + (x+1)%lx;          //(x+1,y  ) sublattice A
            int k = ((y+1)%ly)*lx + (x+1)%lx;          //(x+1,y+1) sublattice A

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 3;
            m->bond2weight[n] = lambda;
            m->sweight += lambda;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            m->bond2index[n*mhnspin+2] = k;
            n++;
        }
    }
    // graph name : reference_vonf
    // range [3*lx*ly, 4*lx*ly)
    // location   : sublattice B
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx          + x + lx*ly;         //(x  ,y  ) sublattice B
            int j = ((y+1)%ly)*lx + x + lx*ly;         //(x  ,y+1) sublattice B
            int k = ((y+1)%ly)*lx + (x+1)%lx + lx*ly;  //(x+1,y+1) sublattice B

            m->bond2type[n]   = 1;
            m->bond2hNspin[n] = 3;
            m->bond2weight[n] = lambda;
            m->sweight += lambda;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            m->bond2index[n*mhnspin+2] = k;
            n++;
        }
    }
    // graph name : gauss_law
    // range [4*lx*ly, 5*lx*ly)
    // location   : sublattice A1
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx        + x;                   //(x  ,y  ) sublattice A
            int j = y*lx + (x+1)%lx;                   //(x+1,y  ) sublattice A

            m->bond2type[n]   = 2;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0;
            m->sweight += 0;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            n++;
        }
    }
    // graph name : gauss_law
    // range [5*lx*ly, 6*lx*ly)
    // location   : sublattice A2
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx          + x;                 //(x  ,y  ) sublattice A
            int j = ((y+1)%ly)*lx + x;                 //(x  ,y+1) sublattice A

            m->bond2type[n]   = 2;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0;
            m->sweight += 0;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            n++;
        }
    }
    // graph name : gauss_law
    // range [6*lx*ly, 7*lx*ly)
    // location   : sublattice A3
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx          + x;                 //(x  ,y  ) sublattice A
            int j = ((y+1)%ly)*lx + (x+1)%lx;          //(x+1,y+1) sublattice A

            m->bond2type[n]   = 2;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0;
            m->sweight += 0;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            n++;
        }
    }
    // graph name : gauss_law
    // range [7*lx*ly, 8*lx*ly)
    // location   : sublattice B1
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx        + x + lx*ly;           //(x  ,y  ) sublattice B
            int j = y*lx + (x+1)%lx + lx*ly;           //(x+1,y  ) sublattice B

            m->bond2type[n]   = 2;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0;
            m->sweight += 0;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            n++;
        }
    }
    // graph name : gauss_law
    // range [8*lx*ly, 9*lx*ly)
    // location   : sublattice B2
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx          + x + lx*ly;         //(x  ,y  ) sublattice B
            int j = ((y+1)%ly)*lx + x + lx*ly;         //(x  ,y+1) sublattice B

            m->bond2type[n]   = 2;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0;
            m->sweight += 0;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            n++;
        }
    }
    // graph name : gauss_law
    // range [9*lx*ly, 10*lx*ly)
    // location   : sublattice B3
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            int i = y*lx          + x + lx*ly;         //(x  ,y  ) sublattice B
            int j = ((y+1)%ly)*lx + (x+1)%lx + lx*ly;  //(x+1,y+1) sublattice B

            m->bond2type[n]   = 2;
            m->bond2hNspin[n] = 2;
            m->bond2weight[n] = 0;
            m->sweight += 0;
            m->bond2index[n*mhnspin+0] = i;
            m->bond2index[n*mhnspin+1] = j;
            n++;
        }
    }

    for(int i=0;i<16;i++) {
        m->link[i] = link_rule_temp_trans[i];
    }
    for(int i=0;i<12;i++) {
        m->link[1*4*mhnspin+i] = link_rule_reference_conf[i];
    }
    for(int i=0;i<8;i++) {
        m->link[2*4*mhnspin+i] = link_rule_gauss_law[i];
    }

    m->insert[0] = insert_rule_temp_trans;
    m->insert[1] = insert_rule_reference_conf;
    m->insert[2] = insert_rule_gauss_law;

    if(n!=nbond) {
        printf("quantum_link_triangular_lattice : something went wrong!\n");
    }

    m->nsite = nsite;
    m->nbond = nbond;
    m->mhnspin = mhnspin;
    create_cmf(m->cmf,m->bond2weight,nbond);

    return m;
}

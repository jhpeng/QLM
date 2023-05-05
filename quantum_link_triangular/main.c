#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "model.h"
#include "update.h"
#include "stats.h"


int uniform_sequence_sampling(double lam, double* sequence, int scap, gsl_rng* rng) {
    int length=0;
    double k = -log(gsl_rng_uniform_pos(rng))/lam;
    while((k<1.0) && (length<scap)) {
        sequence[length] = k;

        k -= log(gsl_rng_uniform_pos(rng))/lam;
        length++;
    }

    return length;
}

static int gausslaw_placeholder_cap=0;
static double* gausslaw_tau_placeholder=NULL;
static int* gausslaw_bond_placeholder=NULL;
void gauss_law_sublattice_A(world_line* w, model* m, int lx, int ly) {
    int sxy = lx*ly;

    if(gausslaw_placeholder_cap==0) {
        gausslaw_placeholder_cap = sxy;
        gausslaw_tau_placeholder = (double*)malloc(sizeof(double)*gausslaw_placeholder_cap);
        gausslaw_bond_placeholder = (int*)malloc(sizeof(int)*gausslaw_placeholder_cap);
    }

    int bond;
    int n=0;
    int s[3];
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            s[0] = (w->istate)[y*lx + x + sxy];                   // (x  , y  ) for sublattice B
            s[1] = (w->istate)[y*lx + (x+1)%lx + sxy];            // (x+1, y  ) for sublattice B
            s[2] = (w->istate)[((y+1)%ly)*lx + (x+1)%lx + sxy];   // (x+1, y+1) for sublattice B

            if(!((s[0] == s[1]) && (s[1] == s[2]))) {
                if(s[0] == s[1]) {
                    bond = ((y+1)%ly)*lx + x + 4*sxy;             // (x  , y+1) for bond A1
                } else if(s[1] == s[2]) {
                    bond = y*lx + x + 5*sxy;                      // (x  , y  ) for bond A2
                } else {
                    bond = y*lx + x + 6*sxy;                      // (x  , y  ) for bond A3
                }

                gausslaw_bond_placeholder[n] = bond;
                gausslaw_tau_placeholder[n] = 0;
                n++;
            }
        }
    }

    insert_vertices(w, m, gausslaw_tau_placeholder, gausslaw_bond_placeholder, n);
}

void gauss_law_sublattice_B(world_line* w, model* m, int lx, int ly) {
    int sxy = lx*ly;

    if(gausslaw_placeholder_cap==0) {
        gausslaw_placeholder_cap = sxy;
        gausslaw_tau_placeholder = (double*)malloc(sizeof(double)*gausslaw_placeholder_cap);
        gausslaw_bond_placeholder = (int*)malloc(sizeof(int)*gausslaw_placeholder_cap);
    }

    int bond;
    int n=0;
    int s[3];
    for(int y=0;y<ly;y++) {
        for(int x=0;x<lx;x++) {
            s[0] = (w->istate)[y*lx + x];                         // (x  , y  ) for sublattice A
            s[1] = (w->istate)[((y+1)%ly)*lx + x];                // (x  , y+1) for sublattice A
            s[2] = (w->istate)[((y+1)%ly)*lx + (x+1)%lx];         // (x+1, y+1) for sublattice A

            if(!((s[0] == s[1]) && (s[1] == s[2]))) {
                if(s[1] == s[2]) {
                    bond = y*lx + x + 7*sxy;                      // (x  , y+1) for bond B1
                } else if(s[0] == s[1]) {
                    bond = y*lx + (x+1)%lx + 8*sxy;               // (x  , y  ) for bond B2
                } else {
                    bond = y*lx + x + 9*sxy;                      // (x  , y  ) for bond B3
                }

                gausslaw_bond_placeholder[n] = bond;
                gausslaw_tau_placeholder[n] = 0;
                n++;
            }
        }
    }

    insert_vertices(w, m, gausslaw_tau_placeholder, gausslaw_bond_placeholder, n);
}

static int sequence_placeholder_cap=0;
static double* sequence_tau_placeholder=NULL;
static int* sequence_bond_placeholder=NULL;
void update_sublattice_A(world_line* w, model* m, int lx, int ly, double lambda, gsl_rng* rng) {
    int sxy = lx*ly;
    double density = (w->beta)*(1+lambda)*(m->nsite)*0.5;

    if(sequence_placeholder_cap==0) {
        sequence_placeholder_cap = 10*((int)density);
        sequence_tau_placeholder = (double*)malloc(sizeof(double)*sequence_placeholder_cap);
        sequence_bond_placeholder = (int*)malloc(sizeof(int)*sequence_placeholder_cap);
    }

    int length = uniform_sequence_sampling(density, sequence_tau_placeholder, sequence_placeholder_cap, rng);
    while(length==sequence_placeholder_cap) {
        sequence_placeholder_cap = 2*sequence_placeholder_cap;

        free(sequence_tau_placeholder);
        free(sequence_bond_placeholder);
        sequence_tau_placeholder = (double*)malloc(sizeof(double)*sequence_placeholder_cap);
        sequence_bond_placeholder = (int*)malloc(sizeof(int)*sequence_placeholder_cap);

        length = uniform_sequence_sampling(density, sequence_tau_placeholder, sequence_placeholder_cap, rng);
    }

    for(int i=0;i<length;i++) {
        if(gsl_rng_uniform_pos(rng)*(1+lambda)<1) {
            sequence_bond_placeholder[i] = (int)(gsl_rng_uniform_pos(rng)*sxy) + sxy;
        } else {
            sequence_bond_placeholder[i] = (int)(gsl_rng_uniform_pos(rng)*sxy) + 2*sxy;
        }
    }

    // removing vetices without kink
    remove_vertices(w);

    // insert vertices on sublattice A from given proposal
    insert_vertices(w, m, sequence_tau_placeholder, sequence_bond_placeholder, length);

    // apply gauss law on sublattice A
    gauss_law_sublattice_A(w, m, lx, ly);

    // clustering and flip the cluster on sublattice A
    clustering(w, m);
    flip_cluster_sublattice_A(w, m, rng);
}

void update_sublattice_B(world_line* w, model* m, int lx, int ly, double lambda, gsl_rng* rng) {
    int sxy = lx*ly;
    double density = (w->beta)*(1+lambda)*(m->nsite)*0.5;

    if(sequence_placeholder_cap==0) {
        sequence_placeholder_cap = 10*((int)density);
        sequence_tau_placeholder = (double*)malloc(sizeof(double)*sequence_placeholder_cap);
        sequence_bond_placeholder = (int*)malloc(sizeof(int)*sequence_placeholder_cap);
    }

    int length = uniform_sequence_sampling(density, sequence_tau_placeholder, sequence_placeholder_cap, rng);
    while(length==sequence_placeholder_cap) {
        sequence_placeholder_cap = 2*sequence_placeholder_cap;

        free(sequence_tau_placeholder);
        free(sequence_bond_placeholder);
        sequence_tau_placeholder = (double*)malloc(sizeof(double)*sequence_placeholder_cap);
        sequence_bond_placeholder = (int*)malloc(sizeof(int)*sequence_placeholder_cap);

        length = uniform_sequence_sampling(density, sequence_tau_placeholder, sequence_placeholder_cap, rng);
    }

    for(int i=0;i<length;i++) {
        if(gsl_rng_uniform_pos(rng)*(1+lambda)<1) {
            sequence_bond_placeholder[i] = (int)(gsl_rng_uniform_pos(rng)*sxy);
        } else {
            sequence_bond_placeholder[i] = (int)(gsl_rng_uniform_pos(rng)*sxy) + 3*sxy;
        }
    }

    // removing vetices without kink
    remove_vertices(w);

    // insert vertices on sublattice B from given proposal
    insert_vertices(w, m, sequence_tau_placeholder, sequence_bond_placeholder, length);

    // apply gauss law on sublattice B
    gauss_law_sublattice_B(w, m, lx, ly);

    // clustering and flip the cluster on sublattice B
    clustering(w, m);
    flip_cluster_sublattice_B(w, m, rng);
}

static int highvar2charge(int* h){
    int charge=0;

    if(h[0]==h[1]) charge+=1;
    else charge-=1;
    if(h[2]==h[3]) charge+=1;
    else charge-=1;
    if(h[4]==h[5]) charge+=1;
    else charge-=1;
    if(h[1]==h[2]) charge-=1;
    else charge+=1;
    if(h[3]==h[4]) charge-=1;
    else charge+=1;
    if(h[5]==h[0]) charge-=1;
    else charge+=1;

    return charge;
}

void counting_charge(int* state, int lx, int ly){
    int indices[6];

    int h[6];
    int charge;



    printf("------------ counting charge --------------\n");
    for(int y=0; y<ly; y++) {
        for(int x=0; x<lx; x++) {
            indices[0] = y*lx+x;
            indices[2] = ((y+1)%ly)*lx+x;
            indices[4] = ((y+1)%ly)*lx+(x+1)%lx;

            indices[1] = y*lx+x                 + lx*ly;
            indices[5] = y*lx+(x+1)%lx          + lx*ly;
            indices[3] = ((y+1)%ly)*lx+(x+1)%lx + lx*ly;

            for(int i=0;i<6;i++)
                h[i] = state[indices[i]];

            charge = highvar2charge(h);

            if(charge!=0)
                printf("(x, y)=(%d, %d) state=(%d %d %d %d %d %d) charge=%d\n",x,y,h[0],h[1],h[2],h[3],h[4],h[5],charge);
        }
    }
}

void initial_state_with_charge(int* state, int lx, int ly, int distance){
    int index1,index2;
    int x,y;

    if(distance > lx/2){
        printf("initial_state_with_charge : the distance can not larger than x/2!\n");
        exit(1);
    }

    for(y=0;y<ly;y++){
        for(x=0;x<lx;x++){
            state[y*lx+x]       = 1;
            state[y*lx+x+lx*ly] = 1;
        }
    }

    y = ly/2-1;
    for(x=(lx-distance)/2;x<(lx+distance)/2;x++){
        index1 = ((y+1)%ly)*lx+(x+1)%lx;

        index2 = y*lx+(x+1)%lx + lx*ly;

        state[index1] = 0;
        state[index2] = 0;
    }

    counting_charge(state,lx,ly);
}

double measure_ma(int* state, int lx, int ly) {
    int sxy = lx*ly;

    int ma=0;
    for(int i=0; i<sxy; i++) {
        ma+=state[i];
    }

    ma = abs(ma*2 - sxy);

    return (double)ma;
}

double measure_mb(int* state, int lx, int ly) {
    int sxy = lx*ly;

    int mb=0;
    for(int i=sxy; i<2*sxy; i++) {
        mb+=state[i];
    }

    mb = abs(mb*2 - sxy);

    return (double)mb;
}

int main(int argc, char** argv) {
    int Lx = 8;
    int Ly = 8;
    double Lambda = 1.0;
    double Beta = 10;
    int distance = 0;
    int Nsample=4000;
    int Nthermal=10000;
    unsigned long int Seed = 1;

    if(argc>1){
        Lx = atoi(argv[1]);
        Ly = atoi(argv[2]);
        Lambda = atof(argv[3]);
        Beta = atof(argv[4]);
        distance = atoi(argv[5]);
        Nthermal = atoi(argv[6]);
        Nsample = atoi(argv[7]);
        Seed = atoi(argv[8]);
    }

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,Seed);

    model* m = quantum_link_triangular_lattice(Lx,Ly,Lambda);

    int Nsite = m->nsite;
    int Mhnspin = m->mhnspin;
    world_line* w = malloc_world_line(1000,2*Mhnspin,Nsite);

    w->beta = Beta;

    initial_state_with_charge(w->istate, Lx, Ly, distance);

    for(int i=0;i<Nthermal;i++) {
        update_sublattice_A(w, m, Lx, Ly, Lambda, rng);
        update_sublattice_B(w, m, Lx, Ly, Lambda, rng);

        counting_charge(w->istate, Lx, Ly);
    }

    double ma1, mb1, ma2, mb2, ma4, mb4;
    estimator* est_ma1 = malloc_estimator(1024,"ma1");
    estimator* est_ma2 = malloc_estimator(1024,"ma2");
    estimator* est_ma4 = malloc_estimator(1024,"ma4");
    estimator* est_mb1 = malloc_estimator(1024,"mb1");
    estimator* est_mb2 = malloc_estimator(1024,"mb2");
    estimator* est_mb4 = malloc_estimator(1024,"mb4");

    int check_point=0;
    for(int i=0; i<Nsample; i++) {
        update_sublattice_A(w, m, Lx, Ly, Lambda, rng);
        update_sublattice_B(w, m, Lx, Ly, Lambda, rng);

        ma1 = measure_ma(w->istate, Lx, Ly);
        mb1 = measure_mb(w->istate, Lx, Ly);
        ma2 = ma1*ma1;
        mb2 = mb1*mb1;
        ma4 = ma2*ma2;
        mb4 = mb2*mb2;

        append_estimator(est_ma1,ma1);
        append_estimator(est_ma2,ma2);
        append_estimator(est_ma4,ma4);
        append_estimator(est_mb1,mb1);
        append_estimator(est_mb2,mb2);
        append_estimator(est_mb4,mb4);

        // save checkpoint
        check_point++;
        if(check_point>=1000) {
            print_detail(est_ma1);
            print_detail(est_mb1);

            save_estimator(est_ma1);
            save_estimator(est_ma2);
            save_estimator(est_ma4);
            save_estimator(est_mb1);
            save_estimator(est_mb2);
            save_estimator(est_mb4);

            check_point=0;
        }
    }

    save_estimator(est_ma1);
    save_estimator(est_ma2);
    save_estimator(est_ma4);
    save_estimator(est_mb1);
    save_estimator(est_mb2);
    save_estimator(est_mb4);

    free_estimator(est_ma1);
    free_estimator(est_ma2);
    free_estimator(est_ma4);
    free_estimator(est_mb1);
    free_estimator(est_mb2);
    free_estimator(est_mb4);
    free_world_line(w);
    free_model(m);

    return 0;
}

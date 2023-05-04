#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "model.h"
#include "update.h"


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

    insert_vertices(w, m, sequence_tau_placeholder, sequence_bond_placeholder, length);

    // apply gauss law on sublattice B
    gauss_law_sublattice_B(w, m, lx, ly);

    // clustering and flip the cluster on sublattice B
    clustering(w, m);
    flip_cluster_sublattice_B(w, m, rng);
}


void main() {
    int Lx = 4;
    int Ly = 4;
    double Beta = 1.2;
    double Lambda = 0.6;
    int Nblock=1000;
    int Nsample=1000;
    int Nthermal=100000;
    unsigned long int seed = 3984793;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    model* m = quantum_link_triangular_lattice(Lx,Ly,lambda);

    int Nsite = m->nsite;
    int Mhnspin = m->mhnspin;
    world_line* w = malloc_world_line(100000,2*Mhnspin,Nsite);

    w->beta = Beta;
}

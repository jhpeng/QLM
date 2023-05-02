#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "dtype.h"
#include "union_find.h"

void remove_vertices(world_line* w) {
    vertex* v;

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    int checindex_old_delete;
    int i,j,index_old;
    index_old=0;
    for(i=0;i<w->nvertices;i++) {
        v = &(sequence1[i]);
        checindex_old_delete = 1;

        for(j=0;j<(v->hNspin);j++) {
            if(v->state[j]!=v->state[j+v->hNspin])
                checindex_old_delete = 0;
        }

        if(!checindex_old_delete) {
            copy_vertex(&(sequence2[index_old]),&(sequence1[i]));
            index_old++;
        }
    }
    w->nvertices = index_old;
    w->flag = !(w->flag);
}

void insert_vertices_sublattice(world_line* w, model* m, double* insert_taus, int* insert_bonds, int insert_len) {
    double sum_weight_sublattice = niste*0.5*(1+lambda);
    double lam = sum_weight_sublattice_A*(w->beta);

    int length = insert_len+w->nvertices;
    realloc_world_line(w,length);

    vertex* sequence1 = w->sequenceB;
    vertex* sequence2 = w->sequenceA;
    if(w->flag) {
        sequence1 = w->sequenceA;
        sequence2 = w->sequenceB;
    }

    int* pstate = w->pstate;
    int  nsite = w->nsite;

    for(int i=0;i<nsite;i++) pstate[i] = w->istate[i];

    vertex* v;
    int n,i,index_old,i_site,index;
    double tau_run,tau_propose;

    int mhnspin = m->mhnspin;
    int lstate[mhnspin];

    index_old=0;
    n=0;

    tau_run = 0;
    if(w->nvertices!=0) tau_run = (sequence1[0]).tau;

    for(i=0;i<insert_len;i++) {
        tau_propose = insert_taus[i];

        while((tau_run<tau_propose) && (index_old<(w->nvertices))) {
            v = &(sequence1[index_old]);
            for(i_site=0;i_site<(v->hNspin);i_site++) {
                index = m->bond2index[v->bond*mhnspin+i_site];
                pstate[index] = v->state[v->hNspin+i_site];
            }

            copy_vertex(&(sequence2[n]),v);
            n++;
            index_old++;

            if(index_old<(w->nvertices)){
                tau_run = (sequence1[index_old]).tau;
            }
        }

        int bond     = insert_bonds[i];
        int t        = m->bond2type[bond];
        int hNspin   = m->bond2hNspin[bond];
        insert_rule rule = m->insert[t];

        for(i_site=0;i_site<hNspin;i_site++) { 
            index = m->bond2index[bond*mhnspin+i_site];
            lstate[i_site] = pstate[index];
        }

        if(rule(lstate)) {
            (sequence2[n]).tau    = tau_propose;
            (sequence2[n]).bond   = bond;
            (sequence2[n]).hNspin = hNspin;

            for(i_site=0;i_site<hNspin;i_site++) {
                (sequence2[n]).state[i_site]        = lstate[i_site];
                (sequence2[n]).state[i_site+hNspin] = lstate[i_site];
            }

            n++;
        }
    }

    while(index_old<(w->nvertices)) {
        copy_vertex(&(sequence2[n]),&(sequence1[index_old]));
        n++;
        index_old++;
    }

    w->nvertices = n;
    w->flag = !(w->flag);
}

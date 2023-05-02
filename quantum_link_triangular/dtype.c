#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dtype.h"

model* malloc_model(int nsite, int nbond, int mhnspin) {
    model* m = (model*)malloc(sizeof(model));

    m->bond2type   = (int*)malloc(sizeof(int)*nbond);
    m->bond2hNspin = (int*)malloc(sizeof(int)*nbond);
    m->bond2weight = (double*)malloc(sizeof(double)*nbond);
    m->bond2index  = (int*)malloc(sizeof(int)*nbond*mhnspin);
    m->link        = (int*)malloc(sizeof(int)*mhnspin*4*20);
    m->insert      = (insert_rule*)malloc(sizeof(insert_rule)*20);
    m->cmf         = (double*)malloc(sizeof(double)*nbond);

    m->nsite = nsite;
    m->nbond = nbond;
    m->mhnspin = mhnspin;
    m->sweight = 0;

    if(1) {
        printf("-------------------------------------------\n");
        printf("#\tmemory allocate : model\n");
        printf("# nsite : %d | nbond : %d | mhnspin : %d\n",nsite,nbond,mhnspin);
        printf("# bond2type   : %zu bytes\n",sizeof(int)*nbond);
        printf("# bond2hNspin : %zu bytes\n",sizeof(int)*nbond);
        printf("# bond2weight : %zu bytes\n",sizeof(double)*nbond);
        printf("# bond2index  : %zu bytes\n",sizeof(int)*nbond*mhnspin);
        printf("# link        : %zu bytes\n",sizeof(int)*mhnspin*4*20);
        printf("# insert      : %zu bytes\n",sizeof(insert_rule)*20);
        printf("# cmf         : %zu bytes\n",sizeof(double)*nbond);
        printf("-------------------------------------------\n");

    }

    return m;
}

void free_model(model* m) {
    free(m->bond2type);
    free(m->bond2hNspin);
    free(m->bond2weight);
    free(m->bond2index);
    free(m->link);
    free(m->insert);
    free(m->cmf);
    free(m);
}

void copy_vertex(vertex* dist, vertex* src) {
    dist->tau = src->tau;
    dist->bond = src->bond;
    dist->hNspin = src->hNspin;
    for(int i=0;i<2*src->hNspin;i++) {
        dist->state[i] = src->state[i];
    }
}

world_line* malloc_world_line(int length, int mnspin, int nsite) {
    world_line* w = (world_line*)malloc(sizeof(world_line));

    w->sequenceA = (vertex*)malloc(sizeof(vertex)*length);
    w->sequenceB = (vertex*)malloc(sizeof(vertex)*length);
    w->cluster  = (int*)malloc(sizeof(int)*length*mnspin);
    w->weight   = (int*)malloc(sizeof(int)*length*mnspin);
    w->istate    = (int*)malloc(sizeof(int)*nsite);
    w->pstate    = (int*)malloc(sizeof(int)*nsite);
    w->first    = (int*)malloc(sizeof(int)*nsite);
    w->last     = (int*)malloc(sizeof(int)*nsite);

    w->length = length;
    w->nvertices = 0;
    w->mnspin = mnspin;
    w->flag = 0;
    w->nsite = nsite;

    return w;
}

void free_world_line(world_line* w) {
    free(w->sequenceA);
    free(w->sequenceB);
    free(w->cluster);
    free(w->weight);
    free(w->istate);
    free(w->pstate);
    free(w->first);
    free(w->last);
    free(w);
}

void realloc_world_line(world_line* w, int length) {
    if(length > (w->length)) {
        vertex* sequence;
        sequence = (vertex*)malloc(sizeof(vertex)*length);
        for(int i=0;i<w->nvertices;i++) 
            copy_vertex(&(sequence[i]),&(w->sequenceA[i]));

        free(w->sequenceA);
        w->sequenceA = sequence;

        sequence = (vertex*)malloc(sizeof(vertex)*length);
        for(int i=0;i<w->nvertices;i++) 
            copy_vertex(&(sequence[i]),&(w->sequenceB[i]));

        free(w->sequenceB);
        w->sequenceB = sequence;

        free(w->cluster);
        free(w->weight);
        w->cluster = (int*)malloc(sizeof(int)*length*(w->mnspin));
        w->weight  = (int*)malloc(sizeof(int)*length*(w->mnspin));

        w->length = length;
    }
}

world_line_omp* malloc_world_line_omp(int cap, int mnspin, int nsite, int nthread) {
    world_line_omp* w = (world_line_omp*)malloc(sizeof(world_line_omp));
    
    w->sequenceA = (vertex**)malloc(sizeof(vertex*)*nthread);
    w->sequenceB = (vertex**)malloc(sizeof(vertex*)*nthread);
    for(int i=0;i<nthread;i++) {
        w->sequenceA[i] = (vertex*)malloc(sizeof(vertex)*cap);
        w->sequenceB[i] = (vertex*)malloc(sizeof(vertex)*cap);
    }

    w->flag = (int*)malloc(sizeof(int)*nthread);
    w->cap  = (int*)malloc(sizeof(int)*nthread);
    w->len  = (int*)malloc(sizeof(int)*nthread);
    for(int i=0;i<nthread;i++) {
        w->flag[i] = 1;
        w->cap[i]  = cap;
        w->len[i]  = 0;
    }

    w->cluster = (int*)malloc(sizeof(int)*cap*mnspin*nthread);
    w->weight  = (int*)malloc(sizeof(int)*cap*mnspin*nthread);

    for(int i=0;i<cap*mnspin*nthread;i++) {
        w->cluster[i] = -1;
        w->weight[i]  = -1;
    }

    w->istate = (int*)malloc(sizeof(int)*nsite*nthread);
    w->pstate = (int*)malloc(sizeof(int)*nsite*nthread);
    w->first  = (int*)malloc(sizeof(int)*nsite*nthread);
    w->last   = (int*)malloc(sizeof(int)*nsite*nthread);

    w->insert_seq = (double**)malloc(sizeof(double*)*nthread);
    w->insert_bond = (int**)malloc(sizeof(int*)*nthread);
    w->insert_cap = (int*)malloc(sizeof(int)*nthread);
    w->insert_len = (int*)malloc(sizeof(int)*nthread);
    for(int i=0;i<nthread;i++) {
        w->insert_seq[i] = (double*)malloc(sizeof(double)*cap);
        w->insert_bond[i] = (int*)malloc(sizeof(int)*cap);
        w->insert_cap[i] = cap;
        w->insert_len[i] = 0;
    }

    w->nthread = nthread;
    w->mcap = cap;
    w->csize = cap*mnspin*nthread;
    w->mnspin = mnspin;
    w->nsite = nsite;

    if(1) {
        printf("-------------------------------------------\n");
        printf("#\tmemory allocate : world_line_omp\n");
        printf("# nsite : %d | mcap : %d | mnspin : %d | nthread %d\n",nsite,cap,mnspin,nthread);
        printf("# sequenceA   : %zu bytes\n",sizeof(vertex)*nthread*cap);
        printf("# sequenceB   : %zu bytes\n",sizeof(vertex)*nthread*cap);
        printf("# cluster     : %zu bytes\n",sizeof(int)*cap*mnspin*nthread);
        printf("# weight      : %zu bytes\n",sizeof(int)*cap*mnspin*nthread);
        printf("# istate      : %zu bytes\n",sizeof(int)*nsite*nthread);
        printf("# pstate      : %zu bytes\n",sizeof(int)*nsite*nthread);
        printf("# first       : %zu bytes\n",sizeof(int)*nsite*nthread);
        printf("# last        : %zu bytes\n",sizeof(int)*nsite*nthread);
        printf("# insert_seq  : %zu bytes\n",sizeof(double)*nthread*cap);
        printf("# insert_bond : %zu bytes\n",sizeof(int)*nthread*cap);
        printf("-------------------------------------------\n");

    }

    return w;
}

void free_world_line_omp(world_line_omp* w) {
    int nthread = w->nthread;

    for(int i=0;i<nthread;i++) {
        free(w->sequenceA[i]);
        free(w->sequenceB[i]);
        free(w->insert_seq[i]);
        free(w->insert_bond[i]);
    }

    free(w->sequenceA);
    free(w->sequenceB);
    free(w->flag);
    free(w->cap);
    free(w->len);
    free(w->cluster);
    free(w->weight);
    free(w->istate);
    free(w->pstate);
    free(w->first);
    free(w->last);
    free(w->insert_seq);
    free(w->insert_bond);
    free(w->insert_cap);
    free(w->insert_len);

    free(w);
}

void realloc_world_line_omp_vertex(world_line_omp* w, int cap, int i_thread) {
    if(w->cap[i_thread]<cap) {
        vertex* sequenceA = (vertex*)malloc(sizeof(vertex)*cap);
        vertex* sequenceB = (vertex*)malloc(sizeof(vertex)*cap);
        vertex* old_seqA = w->sequenceA[i_thread];
        vertex* old_seqB = w->sequenceB[i_thread];

        for(int i=0;i<(w->len[i_thread]);i++) {
            copy_vertex(&(sequenceA[i]),&(old_seqA[i]));
            copy_vertex(&(sequenceB[i]),&(old_seqB[i]));
        }

        free(w->sequenceA[i_thread]);
        free(w->sequenceB[i_thread]);
        w->sequenceA[i_thread] = sequenceA;
        w->sequenceB[i_thread] = sequenceB;

        w->cap[i_thread] = cap;

        if(w->mcap<cap) w->mcap = cap;
    }
}

void realloc_world_line_omp_cluster(world_line_omp* w) {
    if((w->csize)<((w->mcap)*(w->mnspin)*(w->nthread))) {
        free(w->cluster);
        free(w->weight);

        w->cluster = (int*)malloc(sizeof(int)*(w->mcap)*(w->mnspin)*(w->nthread));
        w->weight  = (int*)malloc(sizeof(int)*(w->mcap)*(w->mnspin)*(w->nthread));

        w->csize = (w->mcap)*(w->mnspin)*(w->nthread);

        for(int i=0;i<(w->csize);i++) {
            w->cluster[i] = -1;
            w->weight[i]  = -1;
        }
    }
}

void realloc_world_line_omp_insert(world_line_omp* w, int cap, int i_thread) {
    if(w->insert_cap[i_thread]<cap) {
        double* seq = (double*)malloc(sizeof(double)*cap);
        int* bond = (int*)malloc(sizeof(int)*cap);

        for(int i=0;i<(w->insert_len[i_thread]);i++) {
            seq[i]  = (w->insert_seq[i_thread])[i];
            bond[i] = (w->insert_bond[i_thread])[i];
        }

        free(w->insert_seq[i_thread]);
        free(w->insert_bond[i_thread]);

        w->insert_seq[i_thread] = seq;
        w->insert_bond[i_thread] = bond;
        w->insert_cap[i_thread] = cap;
    }
}

estimator* malloc_estimator(int length, char name[128]) {
    estimator* e = (estimator*)malloc(sizeof(estimator));

    e->samples = (double*)malloc(sizeof(double)*length);
    e->blocks  = (double*)malloc(sizeof(double)*8192);
    strcpy(e->name,name);
    e->length = length;
    e->n = 0;
    e->bsize  = 2;
    e->nblock = 0;

    return e;
}

void free_estimator(estimator* e) {
    free(e->samples);
    free(e->blocks);
    free(e);
}

void realloc_estimator(estimator* e, int length) {
    if(length>(e->length)) {
        double* samples = (double*)malloc(sizeof(double)*length);
        for(int i=0;i<(e->n);i++) 
            samples[i] = e->samples[i];

        free(e->samples);
        e->samples = samples;
        e->length = length;
    }
}

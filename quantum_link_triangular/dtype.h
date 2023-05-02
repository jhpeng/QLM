#ifndef dtype_h
#define dtype_h

typedef int (*insert_rule)(int*);

typedef struct model {
    int* bond2type;
    int* bond2hNspin;
    double* bond2weight;
    double* cmf;
    int* bond2index;
    int* link;
    insert_rule* insert;
    int nsite;
    int nbond;
    int mhnspin;
    double sweight;
} model;

typedef struct vertex {
    double tau;
    int bond;
    int hNspin;
    int state[20];
} vertex;

typedef struct world_line {
    vertex* sequenceA;
    vertex* sequenceB;
    int  length;
    int  nvertices;
    int* cluster;
    int* weight;
    int  mnspin;
    int  flag;
    int* istate;
    int* pstate;
    int* last;
    int* first;
    int  nsite;
    double beta;
} world_line;

typedef struct world_line_omp {
    int nthread;
    vertex** sequenceA;
    vertex** sequenceB;
    int* flag;
    int* cap;
    int* len;
    int mcap;
    int* cluster;
    int* weight;
    int csize;
    int mnspin;
    int nsite;
    int* istate;
    int* pstate;
    int* first;
    int* last;
    double** insert_seq;
    int** insert_bond;
    int* insert_cap;
    int* insert_len;
    double beta;
} world_line_omp;

typedef struct estimator {
    char name[128];
    double* samples;
    int length;
    int n;
    double* blocks;
    int bsize;
    int nblock;
} estimator;

model* malloc_model(
            int nsite, 
            int nbond, 
            int mhnspin);

void free_model(model* m);

void copy_vertex(
            vertex* dist, 
            vertex* src);

world_line* malloc_world_line(
            int length, 
            int mnspin, 
            int nsite);

void free_world_line(world_line* w);

void realloc_world_line(
            world_line* w, 
            int length);

world_line_omp* malloc_world_line_omp(
            int cap, 
            int mnspin, 
            int nsite, 
            int nthread);

void free_world_line_omp(world_line_omp* w);

void realloc_world_line_omp_vertex(
            world_line_omp* w, 
            int cap, 
            int i_thread);

void realloc_world_line_omp_cluster(world_line_omp* w);

void realloc_world_line_omp_insert(
            world_line_omp* w, 
            int cap, 
            int i_thread);

estimator* malloc_estimator(
            int length, 
            char name[128]);

void free_estimator(estimator* e);

void realloc_estimator(
            estimator* e, 
            int length);

#endif

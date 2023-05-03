#ifndef update_h
#define update_h
#include <gsl/gsl_rng.h>

#include "dtype.h"

void remove_vertices(world_line* w);

void insert_vertices(world_line* w, model* m, double* insert_taus, int* insert_bonds, int insert_len);

void clustering(world_line* w, model* m);

void flip_cluster_sublattice_A(world_line* w, model* m, gsl_rng* rng);

void flip_cluster_sublattice_B(world_line* w, model* m, gsl_rng* rng);

#endif

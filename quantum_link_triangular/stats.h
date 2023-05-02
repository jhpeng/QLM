#ifndef stats_h
#define stats_h

#include "dtype.h"

void append_estimator(estimator* e, double sample);

void print_detail(estimator* e);

void save_estimator(estimator* e);

#endif

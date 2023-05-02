#ifndef union_find_h
#define union_find_h

int root(int* p, int v);

void merge(int* p, int* w, int va, int vb);

int root_parallel(int* p, int v);

void merge_parallel(int* p, int* w, int va, int vb);

#endif

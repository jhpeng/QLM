#include <stdio.h>
#include <string.h>
#include <math.h>

#include "dtype.h"

void append_estimator(estimator* e, double sample) {
    if((e->n +1)>(e->length)) {
        realloc_estimator(e,(e->length)*2);
    }

    e->samples[e->n] = sample;
    e->n++;

    if(e->n == (e->bsize)) {
        int i;
        double b = 0;

        for(i=0;i<(e->bsize);i++) 
            b += e->samples[i];

        b = b/(e->bsize);
        e->blocks[e->nblock] = b;
        e->nblock++;

        if(e->nblock==8192) {
            for(i=0;i<4096;i++)
                e->blocks[i] = (e->blocks[2*i] + e->blocks[2*i+1])*0.5;

            e->nblock = 4096;
            e->bsize *= 2;
        }

        e->n=0;
    }
}

static double mean_samples(estimator* e) {
    double mean=0;
    for(int i=0;i<(e->n);i++)
        mean += e->samples[i];

    mean = mean/(e->n);

    return mean;
}

static double std_samples(estimator* e) {
    double mean = mean_samples(e);
    double std = 0;
    double div = 0;

    for(int i=0;i<(e->n);i++) {
        div  = e->samples[i] - mean;
        std += div*div;
    }
    std = sqrt(std/(e->n));

    return std;
}

static double mean_blocks(estimator* e) {
    double mean = 0;
    for(int i=0;i<(e->nblock);i++)
        mean += e->blocks[i];

    mean = mean/(e->nblock);

    return mean;
}

static double std_blocks(estimator* e) {
    double mean = mean_blocks(e);
    double std = 0;
    double div = 0;

    for(int i=0;i<(e->nblock);i++) {
        div  = e->blocks[i] - mean;
        std += div*div;
    }
    std = sqrt(std/(e->nblock));

    return std;
}

static double stderror(estimator* e) {
    return std_blocks(e)/sqrt((double)(e->nblock-1));
}

void print_detail(estimator* e) {
    double mean = mean_blocks(e);
    double err  = stderror(e);
    //double std  = std_samples(e);

    //double std2 = std*std;
    //double tau = err*err/std2*((e->bsize)*(e->nblock));
    //tau = (tau-1)*0.5;

    printf("-------------------------------\n");
    printf("Observable : %s \n",e->name);
    printf("Nblock=%d Bsize=%d \n",e->nblock,e->bsize);
    printf("mean      = %.6e \n",mean);
    printf("std error = %.6e \n",err);
    //printf("std div   = %.6e \n",std);
    //printf("cor time  = %.6e \n",tau);
}

void save_estimator(estimator* e) {
    char filename[128];
    strcpy(filename,e->name);
    strcat(filename,".dat");
    FILE* ofile = fopen(filename,"w");

    fprintf(ofile,"%d %d\n",e->nblock,e->bsize);
    for(int i=0;i<(e->nblock);i++) {
        fprintf(ofile,"%.18e \n",e->blocks[i]);
    }
    fclose(ofile);
}

//
//  data_dist.c
//  
//
//  Created by Eric Strobel on 7/27/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "../process_TECprobe_profiles/UNV/calculate_normalization_factor.h"

#include "data_dist.h"

double mean(double * pop, int n)
{
    int i = 0;
    double sum = 0;
    
    for (i = 0; i < n; i++) {
        sum += pop[i];
    }
    
    return sum/((double)i);
}

double median(double * pop, int n)
{
    double * sorted_pop = NULL;
    double mdn = 0.0;
    double upr = 0.0;
    double lwr = 0.0;
    
    int i = 0;
    
    if ((sorted_pop = malloc(n * sizeof(*sorted_pop))) == NULL) {
        printf("median: error - failed to allocate memory for sorted_pop array. aborting...\n");
        abort();
    }
    
    for (i = 0; i < n; i++) {
        sorted_pop[i] = pop[i];
    }
    
    qsort(&sorted_pop[0], n, sizeof(sorted_pop[0]), cmpfnc);
    
    
    if (n % 2) {
        mdn = sorted_pop[n/2];
    } else {
        upr = sorted_pop[n/2];
        lwr = sorted_pop[(n/2)-1];
        mdn = (upr + lwr)/2;
    }
    
    free(sorted_pop);
    
    return mdn;
}

double stddev(double * pop, int n, int typ)
{
    int i = 0;
    
    double u = mean(pop, n);
    double nmrtr = 0.0;

    for (i = 0; i < n; i++) {
        nmrtr += pow(pop[i]-u, 2);
    }
    
    if (typ == POPULATION) {
        return sqrt(nmrtr/n);
    } else if (typ == SAMPLE) {
        return sqrt(nmrtr/(n-1));
    } else {
        printf("stddev: error - type variable must be POPULATION or SAMPLE. aborting...\n");
        abort();
    }
}

double calc_Zscore(double x, double mean, double sd)
{
    return (x-mean)/sd;
}

double normalCDF(double z)
{
    /* printf("z:     %f\n", z);
    printf("sr2:   %f\n", M_SQRT2);
    printf("z/sr2: %f\t%f\n", z/M_SQRT2, -1.452326/1.41421356);
    printf("erf:   %f\t%f\n", erf(z/M_SQRT2), erf(-1.027174));
    */
    return 0.5 * (1 + erf(z/M_SQRT2));
}





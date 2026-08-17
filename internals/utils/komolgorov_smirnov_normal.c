//
//  komolgorov_smirnov_normal.c
//  
//
//  Created by Eric Strobel on 7/29/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "../process_TECprobe_profiles/UNV/calculate_normalization_factor.h"

#include "data_dist.h"

#include "komolgorov_smirnov_normal.h"

double komolgorov_smirnov_normal(double * ipt, int n)
{
    int ix = 0;
    int i  = 0;
    
    double * x = NULL; //sorted input
    double * z = NULL; //Z-scores of input data
    double * nCDF = NULL;
    
    double u = mean(ipt, n);
    double s = stddev(ipt, n, SAMPLE);
    
    printf("u = %f\n", u);
    printf("s = %f\n", s);
    
    double D =  -DBL_MAX;
    double Dp = 0.0;
    double Dn = 0.0;
    double * Dm = NULL;
    
    if ((x = malloc(n * sizeof(*x))) == NULL) {
        printf("Komolgorov_Smirnov_normal: error - failed to allocate memory for sorted input values. aborting...\n");
        abort();
    }
    
    if ((z = malloc(n * sizeof(*z))) == NULL) {
        printf("Komolgorov_Smirnov_normal: error - failed to allocate memory for Z-scores. aborting...\n");
        abort();
    }
    
    if ((nCDF = malloc(n * sizeof(*nCDF))) == NULL) {
        printf("Komolgorov_Smirnov_normal: error - failed to allocate memory normal CDF values. aborting...\n");
        abort();
    }
    
    for (ix = 0; ix < n; ix++) {
        x[ix] = ipt[ix];
    }
    
    qsort(&x[0], n, sizeof(x[0]), cmpfnc);
    
    /*for (i = -3; i <= 3; i++) {
        printf("cdf %2d: %f\n", i, normalCDF(i));
    }*/
    
    for (ix = 0, i = 1; ix < n; ix++, i++) {
        
        printf("x[%d] = %f\n", ix, x[ix]);
        
        z[ix] = calc_Zscore(x[ix], u, s);
        nCDF[ix] = normalCDF(z[ix]);
        printf("z[%d]    = %f\n", ix, z[ix]);
        printf("nCDF[%d] = %f\n", ix, nCDF[ix]);
        
        Dp = fabs((((double)i)/n) - nCDF[ix]);
        Dn = fabs(nCDF[ix] - ((((double)i)-1)/n));
        
        printf("Dp = %f\nDn = %f\n\n", Dp, Dn);
        
        Dm = (Dp > Dn) ? &Dp : &Dn;
        
        if (*Dm > D) {
            D = *Dm;
        }
        
    }
    
    
    free(x);
    free(z);
    free(nCDF);
    
    return D;
}

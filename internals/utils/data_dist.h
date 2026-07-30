//
//  data_dist.h
//  
//
//  Created by Eric Strobel on 7/27/26.
//

#ifndef data_dist_h
#define data_dist_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define POPULATION 0
#define SAMPLE 1

double mean(double * pop, int n);
double median(double * pop, int n);
double stddev(double * pop, int n, int typ);
double calc_Zscore(double x, double mean, double sd);
double normalCDF(double z);



#endif /* data_dist_h */

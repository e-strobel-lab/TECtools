//
//  shapiro_wilk_test.h
//  
//
//  Created by Eric Strobel on 7/29/26.
//

#ifndef shapiro_wilk_test_h
#define shapiro_wilk_test_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "../process_TECprobe_profiles/UNV/calculate_normalization_factor.h"

#include "data_dist.h"

#define SW_TBL_MN 1
#define SW_p_0p01 1
#define SW_p_0p02 2
#define SW_p_0p05 3
#define SW_p_0p10 4
#define SW_p_0p50 5
#define SW_p_0p90 6
#define SW_p_0p95 7
#define SW_TBL_MX 7

double shapiro_wilk_test(double * pop, int n);
double royston_approx_4to11(double W, double n);
double royston_approx_12to2K(double W, double n);

#endif /* shapiro_wilk_test_h */

//
//  komolgorov_smirnov_normal.h
//  
//
//  Created by Eric Strobel on 7/29/26.
//

#ifndef komolgorov_smirnov_normal_h
#define komolgorov_smirnov_normal_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "../process_TECprobe_profiles/UNV/calculate_normalization_factor.h"

#include "data_dist.h"

double komolgorov_smirnov_normal(double * ipt, int n);

#endif /* komolgorov_smirnov_normal_h */

//
//  calculate_normalization_factor.h
//  
//
//  Created by Eric Strobel on 2/13/24.
//

#ifndef calculate_normalization_factor_h
#define calculate_normalization_factor_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./store_SM2_profile.h"

#define PCT90TH 0
#define PCT95TH 1

/* calcualate_normalization_factor - perform reactivity profile normalization. the normalization procedure used is adapted from shapemapper2 (https://github.com/Weeks-UNC/shapemapper2) */
double calculate_normalization_factor(double * reactivity_list, int len, int pct2use);

/* cmpfunc: simple float comparison function used to sort reactivity array. */
int cmpfnc(const void * a, const void * b);

/* quantile: C implementation of the Python quantile function described here: https://web.archive.org/web/20150906230750/http://adorio-research.org/wordpress/?p=125
   which is used by shapemapper2 during reactivity normalization and is an implementation of
   the Quantile function in R.
 
   From the source webpage: "given an array x, and a quantile value from 0 to 1.0 (q), the
   function will return a value of x which may not be an element of the array such that
   P(X <= x_q) < q, i.e. q * 100 percent of the data are less than or equal to x_q."
 
   The function will abort the program if an unsorted array is provided. I chose this outcome
   rather than attempting to sort the array to avoid complications from arrays that contain NaNs
 */
double quantile(double * x, int n, double q, int qtype);

/* isHQnuc: assess whether reactivity value is high quality */
int isHQnuc(SM2_profile * prf, int nt_ix, int min_depth, double max_bkg);

#endif /* calculate_normalization_factor_h */

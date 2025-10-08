//
//  calculate_normalization_factor.c
//  
//
//  Created by Eric Strobel on 2/13/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"
#include "../../utils/gen_utils.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./store_SM2_profile.h"

#include "calculate_normalization_factor.h"

/* calcualate_normalization_factor - perform reactivity profile normalization. the normalization procedure used is adapted from shapemapper2 (https://github.com/Weeks-UNC/shapemapper2) */
double calculate_normalization_factor(double * reactivity_list, int len, int pct2use)
{
    if (len < 10) { //return 0 if too few reactivities to perform normalization
        printf("calculate_normalization_factor: error - sequence contains too few nucleotides with quality reactivity information for effective normalization factor calculation.\n");
        return 0;
    }
    
    int i = 0; //general purpose index
    
    //IQR threshold variables
    double q25 = 0;
    double q75 = 0;
    double q_limit = 0;
    
    //90th/95th percentile threshold variables
    int pct10 = 0;         //10% of the data
    int pct05 = 0;         //5% of the data
    double limit90 = 0;    //90th percentile limit
    double limit95 = 0;    //95th percentile limit
    double pct_limit = 0;  //percentile threshold reactivity value
    
    double limit = 0;      //limit to use for boxplot normalization
    
    int upr_thresh_ix = 0; //index of the upper threshold limit
    int lwr_thresh_ix = 0; //incex of the lower threshold limit
    
    //sort reactivity list
    qsort(&reactivity_list[0], len, sizeof(reactivity_list[0]), cmpfnc);
    
    //calculate the 90th/95th percentile index and limit
    pct10 = (int)(floor)(len/10);
    pct05 = (int)(floor)(len/20);
    limit90 = reactivity_list[len - 1 - pct10];
    limit95 = reactivity_list[len - 1 - pct05];
    
    //set the percentile limit to use
    if (pct2use == PCT90TH) {
        pct_limit = limit90;
    } else if (pct2use == PCT95TH) {
        pct_limit = limit95;
    } else {
        printf("calculate_normalization_factor: error - unexpected value for pct2use. aborting...");
        abort();
    }
    
    //calculate q_limit
    q25 = quantile(reactivity_list, len, 0.25, 7);
    q75 = quantile(reactivity_list, len, 0.75, 7);
    q_limit = 1.5 * fabs(q25-q75);
    
    //select the limit that removes the least data
    limit = (pct_limit > q_limit) ? pct_limit : q_limit;
    
    //iterate through the reactivity list until thresh2use is exceeded
    //to find the upper threshold index
    for (i = 0; i < len && reactivity_list[i] < limit; i++) {;}
    
    //set upper threshold index
    upr_thresh_ix = i;
    
    //set lower threshold index as the 90th percentile
    //reactivity relative to the upper threshold index
    lwr_thresh_ix = upr_thresh_ix - pct10;
    
    if (lwr_thresh_ix < 0 || lwr_thresh_ix >= len) {
        printf("calculate_normalization_factor: error - lower threshold index (%d) is outside of array bounds (%d, %d). aborting...\n", lwr_thresh_ix, 0, len);
        abort();
    } else if (upr_thresh_ix < 0 || upr_thresh_ix >= len) {
        printf("calculate_normalization_factor: error - upper threshold index (%d) is outside of array bounds (%d, %d). aborting...\n", upr_thresh_ix, 0, len);
        abort();
    }
    
    //calculate the normalization factore by finding the
    //average of the top 10% of reactivities relative to
    //the upper threshold limit (the average of all
    //reactivity values between the upper and lower
    //threshold limits
    
    double rct_sum = 0; //sum of reactivity values, used for averaging
    int denom = 0;      //number of reactivity values between the upper and lower thresholds
    
    for (i = lwr_thresh_ix; i < upr_thresh_ix; i++) { //for values between upr and low thresh
        rct_sum += reactivity_list[i];                //add to reactivity sum
        denom++;                                      //increment number of reactivity vals summed
    }

    return (rct_sum/((double)denom)); //calculate the average
    
    /* test data from original implementation of the quantile function
    double test_array[9] = {11.4, 17.3, 21.3, 25.9, 40.1, 50.5, 60.0, 70.0, 75};
    
    for (i = 1; i <= 9; i++) {
        printf("%d\t%f\n", i, quantile(test_array, 9, 0.35, i, 1));
    }
     expected output:
     1  25.9
     2  25.9
     3  21.3
     4  21.99
     5  24.29
     6  23.6
     7  24.98
     8  24.06
     9  24.1175
     */
}

/* cmpfunc: simple float comparison function used to sort reactivity array. */
int cmpfnc(const void * a, const void * b)
{
    if (*(double*)a - *(double*)b > 0) {
        return 1;
        
    } else if (*(double*)a - *(double*)b < 0) {
        return -1;
        
    } else {
        return 0;
    }
    
}

/* quantile: C implementation of the Python quantile function described here: https://web.archive.org/web/20150906230750/http://adorio-research.org/wordpress/?p=125
   which is used by shapemapper2 during reactivity normalization and is an implementation of
   the Quantile function in R.
 
   From the source webpage: "given an array x, and a quantile value from 0 to 1.0 (q), the
   function will return a value of x which may not be an element of the array such that
   P(X <= x_q) < q, i.e. q * 100 percent of the data are less than or equal to x_q."
 
   The function will abort the program if an unsorted array is provided. I chose this outcome
   rather than attempting to sort the array to avoid complications from arrays that contain NaNs
 */
double quantile(double * x, int n, double q, int qtype)
{
    //check that array is sorted, abort if not
    int ix = 0;
    for (ix = 0; ix < n; ix++) {
        if (ix > 0 && x[ix] < x[ix-1]) {
            printf("quantile: error - data is not sorted. aborting...\n");
            abort();
        }
    }
    
    //quantile parameters
    double params[10][4] = {
        {0, 0, 1, 0},    //inverse empirical distrib.function., R type 1
        {0.5, 0, 1, 0},  //similar to type 1, averaged, R type 2
        {0.5, 0, 0, 0},  //nearest order statistic,(SAS) R type 3
        {0, 0, 0, 1},    //California linear interpolation, R type 4
        {0.5, 0, 0, 1},  //hydrologists method, R type 5
        {0, 1, 0, 1},    //mean-based estimate(Weibull method), (SPSS,Minitab), type 6
        {1, -1, 0, 1},   //mode-based method,(S, S-Plus), R type 7
        {(1.0/3), (1.0/3), 0, 1}, //median-unbiased ,  R type 8
        {(3/8.0), 0.25, 0, 1}     //normal-unbiased, R type 9.
    };
    
    //set quantile parameters to the user-spedified type
    double a = params[qtype-1][0];
    double b = params[qtype-1][1];
    double c = params[qtype-1][2];
    double d = params[qtype-1][3];
    
    double g = 0; //decimal part of fractional array index
    double j = 0; //integer part of fractional array index
    int i = 0;    //variable used when forcing integer part of frac array index from double to int
    
    //calculate a fractional array index for quantile q,
    //separating the integer (j) and decimal (g) components
    g = modf(a + (n + b) * q - 1, &j);
    
    if (j < 0) {         //quantile index is less than first index (0)
        return x[0];     //return value of first member of array x
        
    } else if (j >= n) { //quantile index is greater than the maximum array index
        return x[n-1];   //return value of second to last member of array x
    }
    
    i = (int)(floor(j)); //convert integer part of frac array index from double to int
    
    if (g == 0) {        //fractional array index has no decimal part
        return x[i];     //return value of array x at index i
    } else {             //fractional array index has a decimal part
        return x[i] + (x[i+1] - x[i]) * (c + d * g); //return non-array member value at quantile q
    }
}

/* isHQnuc: assess whether reactivity value is high quality. the criteria for
 HQ nucleotides used in this function match those used by shapemapper2 */
int isHQnuc(SM2_profile * prf, int nt_ix, int min_depth, double max_bkg) {
    
    //exclude left most nucleotide
    if (nt_ix == 0) {
        return 0;
    }
    
    //test whether nucleotide is masked as lowercase
    if (!isupper(prf->sequence[nt_ix])) {
        return 0;
    }
    
    //test whether reactivity is a number
    if (isnan(prf->reactivity_profile[nt_ix])) {
        return 0;
    }
    
    //test whether modified effective depth is sufficient
    if (prf->mod_eff_depth[nt_ix] < min_depth) {
        return 0;
    }
    
    //test whether untreated effective depth is sufficient
    if (prf->chnls.unt && prf->unt_eff_depth[nt_ix] < min_depth) {
        return 0;
    }
    
    //test whether denatured effective depth is sufficient
    if (prf->chnls.den && prf->den_eff_depth[nt_ix] < min_depth) {
        return 0;
    }

    //test that background is less than the maximum
    if (prf->unt_rate[nt_ix] > max_bkg) {
        return 0;
    }
    
    return 1;
}

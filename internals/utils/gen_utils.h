//
//  gen_utils.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef gen_utils_h
#define gen_utils_h

#include <stdio.h>
#include <math.h>
#include <string.h>

/* compare_float: test equality for floating point values a and b by
 testing whether the the absolute value of a-b is less than a user-provided
 precision setting. */
int compare_float(double a, double b, double precision);

/* printf2_scrn_n_fl: print a string to the screen and to a file */
void printf2_scrn_n_fl(FILE * out_fp, char * out_str);

/* chk_filepath_frmt: ensure that filepath ends in '/' char */
void chk_filepath_frmt(char * path);

#endif /* gen_utils_h */

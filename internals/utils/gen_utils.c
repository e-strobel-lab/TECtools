//
//  gen_utils.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "gen_utils.h"


/* compare_float: test equality for floating point values a and b by
 testing whether the the absolute value of a-b is less than a user-provided
 precision setting. */
int compare_float(double a, double b, double precision)
{
    return ((fabs(a-b) < precision) ? 1 : 0);
}


/* printf2_scrn_n_fl: print a string to the screen and to a file */
void printf2_scrn_n_fl(FILE * out_fp, char * out_str)
{
    printf("%s", out_str);
    fprintf(out_fp, "%s",  out_str);
}


/* chk_filepath_frmt: ensure that filepath ends in '/' char */
void chk_filepath_frmt(char * path)
{
    int last_char = 0;
    
    last_char = strlen(path)-1; //set last_char index

    if (path[last_char]  != '/') {
        path[last_char+1] = '/';
        path[last_char+2] = '\0';
    }
}

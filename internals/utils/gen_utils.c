//
//  gen_utils.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

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

/* check_int_str: check that string is composed of digits */
int check_int_str(char * str, int action)
{
    int i = 0;
    
    int valid = 1;
    
    if (!strcmp(str, "nan")) {
        valid = 1;
    } else {
        for (i = 0; str[i] && valid; i++) {
            if (!isdigit(str[i])) {
                valid = 0;
            }
        }
    }
    
    if (action == RETURN_OUTCOME || valid) {
        return valid;
        
    } else if (action == ABORT_FAILURE && !valid) {
        printf("check_int_str: integer string '%s' contains the non-digit character '%c'. aborting...\n", str, str[i]);
        abort();
        
    } else {
        printf("check_int_str: unexpected action code. aborting...\n");
        abort();
    }
}

/* check_float_str: check that string is composed of valid float chars */
int check_float_str(char * str, int action)
{
    int i = 0;
    
    int valid = 1;
    
    if (!strcmp(str, "nan")) {
        valid = 1;
    } else {
        for (i = 0; str[i] && valid; i++) {
            if (!isdigit(str[i]) &&
                str[i] != '.'    &&
                str[i] != '-'    &&
                str[i] != 'e'    &&
                str[i] != 'E') {
                valid = 0;
            }
        }
    }
    
    if (action == RETURN_OUTCOME || valid) {
        return valid;
        
    } else if (action == ABORT_FAILURE && !valid) {
        printf("check_float_str: float string '%s' contains the invalid character '%c'. aborting...\n", str, str[i]);
        abort();
        
    } else {
        printf("check_float_str: unexpected action code. aborting...\n");
        abort();
    }
}

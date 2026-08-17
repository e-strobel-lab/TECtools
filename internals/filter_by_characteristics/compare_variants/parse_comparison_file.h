//
//  parse_comparison_file.h
//  
//
//  Created by Eric Strobel on 6/29/26.
//

#ifndef parse_comparison_file_h
#define parse_comparison_file_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"

#define MAX_CMP_COL_NM  64
#define MAX_COMPARISONS 8

#define SIM 0
#define DIF 1

/* comparison_values: structure for storing values proved by the comparison input file*/
typedef struct comparison_values {
    char * nm;
    int typ;
    int ix;
} comparison_values;

/* parse_comparison file: parse and store comparison details from input comparison values file */
int parse_comparison_file(char * cmp_path, comparison_values ** cmp);

/* count_entries: count number of lines in input comparison file */
int count_entries(char * path, int lpe);

#endif /* parse_comparison_file_h */

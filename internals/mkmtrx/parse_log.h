//
//  parse_log.h
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#ifndef parse_log_h
#define parse_log_h

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"
#include "../global/global_defs.h"

/* parse_log: open and parse log file to get bowtie2 alignment statistics
 for modified and untreated channels */
int parse_log(char * file_loc, alignment_stats * algn);

/* get_leading_digits: store the first encountered
 contiguous string of digits in the digits array */
int get_leading_digits(char * line, char * digits);

/* calc_reads_algned: calculate the total number of reads that were
 assessed and the total number of reads that were aligned */
int calc_reads_algnd(alignment_stats * algn, cotrans_matrix * mtrx);

#endif /* parse_log_h */

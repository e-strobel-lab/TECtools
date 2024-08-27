//
//  assemble_TECprobeLM_data_structs.h
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#ifndef assemble_TECprobeLM_data_structs_h
#define assemble_TECprobeLM_data_structs_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../mkmtrx/cotrans_mtrx.h"
#include "./assemble_TECprobeLM_data_defs.h"
#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"

/* mode_parameters: structure for storing mode-specific variables */
typedef struct mode_parameters {
    int  mod;           //mode code
    char dlm;           //delimiter used in file
    char hdr[MAX_NAME]; //string to search for the target header substring in column headers
    int offset;         //amount to offset array indices when printing nucleotide
} mode_parameters;

/* sample_names: structure for storing sample names */
typedef struct input_data {
    FILE * fp[TOT_SAMPLES][MAX_IPT]; //input file pointers
    char * fn[TOT_SAMPLES][MAX_IPT]; //array to store input file names
    char * sn[TOT_SAMPLES][MAX_IPT]; //array to store input sample names
    configuration_MLT cfg[TOT_SAMPLES][MAX_IPT]; //configs to store parsed sample name attributes
    int cnt[TOT_SAMPLES];            //number of input files for each sample
} input_data;


#endif /* assemble_TECprobeLM_data_structs_h */

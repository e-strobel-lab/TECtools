//
//  parse_VL_sample_name.h
//  
//
//  Created by Eric Strobel on 2/5/24.
//

#ifndef parse_VL_sample_name_h
#define parse_VL_sample_name_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../utils/io_management.h"

#include "../../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"
#include "../../cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_nm.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

/* parsed_sample_name: pointers to sample name sub-strings generated when parsing the sample_name string */
typedef struct parsed_sample_name {
    char *anchr;             //pointer to initial anchor ('_CoTxn_' or '_Equil_') used to parse sample name
    char *rna;               //pointer to rna name string
    char *fld;               //pointer to folding condition string
    char *prb;               //pointer to probe string
    char *lig;               //pointer to ligand string
    char *conc;              //pointer to ligand concentration string
    char *runID[MAX_RUNS];   //pointer to run ID string
    char *run_cnt_S;         //pointer to run count string
    char *SM;                //pointer to smoothing flag string
    char *rst;               //pointer to the rest of the sample name (optional fields)
    char *field[MAX_FIELDS]; //pointer to custom value field
    int run_cnt_I;           //integer value of run count string
    int runID_cnt;           //number of runIDs that have been parsed
    int field_cnt;           //number of custom value fields that have been parsed
} parsed_sample_name;

/* parse_VL_sample_name: parse sample name for attributes that were specified in the TECprobe analysis config */
void parse_VL_sample_name(char * ipt_nm, configuration_MLT * cfg);

/* remove_out_suffix: remove "_out" suffix from sample name */
void remove_out_suffix (char * snm);

/* remove_simple_uffix: remove a string of characters from the end of a TECprobe sample name */
void remove_simple_suffix(char * sn, char * str2rmv);

/* print_parsed_fields: print parsed fields from config values */
void print_parsed_fields(char * ipt_nm, configuration_MLT * cfg);

#endif /* parse_VL_sample_name_h */

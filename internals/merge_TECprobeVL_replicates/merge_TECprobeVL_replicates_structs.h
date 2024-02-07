//
//  merge_TECprobeVL_replicates_structs.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef merge_TECprobeVL_replicates_structs_h
#define merge_TECprobeVL_replicates_structs_h

#include <stdio.h>
#include <dirent.h>

#include "../global/global_defs.h"
#include "../mkmtrx/mkmtrx_defs.h"
#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"

#include "./merge_TECprobeVL_replicates_defs.h"

/* SM2_analysis_directory: pointers to all relevant directories and files of a TECprobe-ML data set */
typedef struct SM2_analysis_directory {
    char prnt_dir_nm[MAX_NAME+1]; //parent directory name
    DIR * prnt;                   //pointer to parent directory
    DIR * tl[MAX_ROW];            //pointer to transcript length directory
    DIR * out[MAX_ROW];           //pointer to SM2 output directory
    FILE * prf[MAX_ROW];          //pointer to profile file
    int min_tl;                   //minimum transcript length
    int max_tl;                   //maximum transcript length
} SM2_analysis_directory;

/* output_files: pointers and file names for merged output files */
typedef struct output_files {
    char out_dir[MAX_NAME];     //output directory name
    FILE * ofp[MAX_ROW];          //output file pointers
    char ofn[MAX_ROW][MAX_NAME];  //output file names
} output_files;

/* sample_names: structure to manage input sample name parsing and merged sample name construction */
typedef struct sample_names {
    char ipt[MAX_RUNS][MAX_NAME];    //input file sample names
    int cnt;                         //number of input file sample names
    char usr[MAX_NAME];              //user-supplied sample name
    char mrg[MAX_NAME];              //merged sample name
    char * sn2use;                   //pointer to sample name to use for output files
    configuration_MLT cfg[MAX_RUNS]; //stores parsed name info
    configuration_MLT mrgd_cfg;      //stores merged name info for automated sample name construction
} sample_names;

#endif /* merge_replicate_SM2_out_structs_h */

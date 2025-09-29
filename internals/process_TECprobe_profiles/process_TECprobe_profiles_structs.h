//
//  process_TECprobe_profiles_structs.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef process_TECprobe_profiles_structs_h
#define process_TECprobe_profiles_structs_h

#include <stdio.h>
#include <dirent.h>

#include "../global/global_defs.h"
#include "../mkmtrx/mkmtrx_defs.h"
#include "../cotrans_preprocessor/run_script_gen/UNV/config_struct.h"
#include "./UNV/store_SM2_profile.h"

#include "./process_TECprobe_profiles_defs.h"

/* SM2_analysis_directory: pointers to all relevant directories and files of a TECprobe-ML data set */
typedef struct SM2_analysis_directory {
    char prnt_dir_nm[MAX_NAME+1]; //parent directory name
    char ** loc;                  //relative filepaths of profiles
    SM2_profile * data;           //SM2 profile data
    int * indx;                   //index table for location/profile indices
    channel_tracker chnls;        //struct to track what channels are present
    int sd_cnt;                   //target analysis subdirectory count
    int outs_cnt;                 //number of output directories opened
    int prfs_cnt;                 //number of profiles found
    int trgt_start;               //target RNA start index
    int min_id;                   //minimum id
    int max_id;                   //maximum id
    int len[2];                   //min and max target length
    int trg_rct_cnt;              //target nucleotide reactivity count
    double cnf;                   //calculated normalization factor
} SM2_analysis_directory;

/* sample_names: structure to manage input sample name parsing and merged sample name construction */
typedef struct sample_names {
    char ipt[MAX_RUNS][MAX_NAME];       //input file sample names
    int cnt;                            //number of input file sample names
    char usr[MAX_NAME];                 //user-supplied sample name
    char mrg[MAX_NAME];                 //merged sample name
    char * sn2use;                      //pointer to sample name to use for output files
    tprobe_configuration cfg[MAX_RUNS]; //stores parsed name info
    tprobe_configuration mrgd_cfg;      //stores merged name info for automated sample name construction
} sample_names;

#endif /* process_TECprobe_profiles_structs_h */

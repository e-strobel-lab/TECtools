//
//  mk_MLT_run_script.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/debug.h"
#include "config_MLT_struct.h"
#include "parse_MLT_config.h"
#include "check_MLT_config.h"
#include "mk_MLT_run_nm.h"
#include "print_MLT_SM2_script.h"

#include "mk_MLT_run_script.h"

const char VL_LM_config_header[36] = "TECprobe-VL/TECprobe-LM config file";
const char SL_config_header[24] = "TECprobe-SL config file";
const char MUX_config_header[25] = "TECprobe-MUX config file";

extern int debug; //flag to run debug mode

/* mk_MLT_run_script: manages shapemapper2 run script generation for multilength cotranscriptional RNA
 structure probing experiments */
int mk_MLT_run_script(FILE * fp_config_MLT)
{
    //declare and initialize config_MLT structure
    //some config_MLT variables are initialized to -1
    configuration_MLT config_MLT = {0};
    config_MLT.min_target_length = -1;
    config_MLT.max_target_length = -1;
    config_MLT.concatenated = -1;
    config_MLT.smoothing = -1;
    config_MLT.cotranscriptional = -1;
    if (debug) {check_cnfgMLT_init(&config_MLT);} //in debug mode, check initialization of config struct
    
    int mode = -1;                    //mode variable that will be set based on config header
    char sample_name[MAX_LINE] = {0}; //name that will be used for shapemapper2 analysis
    
    parse_MLT_config(fp_config_MLT, &config_MLT, &mode); //parse input config file
    check_MLT_config(&config_MLT, mode); //check global variables for each setting to confirm that they were set correctly
    mk_MLT_run_nm(sample_name, &config_MLT);		//construct sample name for shapemapper2
    print_MLT_SM2_script(sample_name, &config_MLT); //generate shapemapper2 run script
    
    printf("%s\n\n", sample_name);    //print final sample name
    
    return 1;
}

/* check initialization: check initialization of config struct */
void check_cnfgMLT_init(configuration_MLT * config_MLT)
{
    printf("\ninitial config values\n");
    printf("shapemapper2_path\t%2d\n",      config_MLT->shpmppr2_path[0]);
    
    printf("input_file_location\t%2d\n",    config_MLT->ipt_file_loc[0]);
    printf("output_file_location\t%2d\n",   config_MLT->out_file_loc[0]);
    
    printf("untreated_read1_prefix\t%2d\n", config_MLT->untreated_read1_prefix[0]);
    printf("untreated_read2_prefix\t%2d\n", config_MLT->untreated_read2_prefix[0]);
    
    printf("modified_read1_prefix\t%2d\n",  config_MLT->modified_read1_prefix[0]);
    printf("modified_read2_prefix\t%2d\n",  config_MLT->modified_read2_prefix[0]);
    
    printf("target_files_location\t%2d\n",  config_MLT->trg_files_loc[0]);
    printf("target_files_prefix\t%2d\n",    config_MLT->trg_files_prfx[0]);
    printf("min_target_length\t%2d\n",      config_MLT->min_target_length);
    printf("max_target_length\t%2d\n",      config_MLT->max_target_length);
    
    printf("input_name\t\t%2d\n",           config_MLT->input_name[0]);
    
    printf("cotranscriptional\t%2d\n",      config_MLT->cotranscriptional);
    printf("chemical_probe\t\t%2d\n",       config_MLT->chemical_probe[0]);
    
    printf("concatenated\t\t%2d\n",         config_MLT->concatenated);
    printf("runID\t\t\t%2d\n",              config_MLT->runID[0][0]);
    printf("smoothing\t\t%2d\n",            config_MLT->smoothing);
    
    printf("ligand_name\t\t%2d\n",          config_MLT->ligand_name[0]);
    printf("ligand_conc\t\t%2d\n",          config_MLT->ligand_conc[0]);
    
    printf("field\t\t\t%2d\n",              config_MLT->field[0][0]);
    printf("run_count\t\t%2d\n",            config_MLT->run_count);
    printf("field_count\t\t%2d\n\n",        config_MLT->field_count);
}

//
//  mk_run_script.c
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
#include "config_struct.h"
#include "parse_config.h"
#include "check_config.h"
#include "../UNV/mk_run_nm.h"
#include "../MLT/print_MLT_SM2_script.h"

#include "mk_run_script.h"

const char VL_LM_config_header[36] = "TECprobe-VL/TECprobe-LM config file";
const char SL_config_header[24] = "TECprobe-SL config file";
const char MUX_config_header[25] = "TECprobe-MUX config file";

extern int debug; //flag to run debug mode

/* mk_run_script: manages shapemapper2 run script generation for multilength cotranscriptional RNA
 structure probing experiments */
int mk_run_script(FILE * fp_config)
{
    //declare and initialize config structure
    //some config variables are initialized to -1
    tprobe_configuration config = {0};
    config.min_target_length = -1;
    config.max_target_length = -1;
    config.concatenated = -1;
    config.smoothing = -1;
    config.cotranscriptional = -1;
    if (debug) {check_cnfgMLT_init(&config);} //in debug mode, check initialization of config struct
    
    int mode = -1;                    //mode variable that will be set based on config header
    char sample_name[MAX_LINE] = {0}; //name that will be used for shapemapper2 analysis
    
    parse_config(fp_config, &config, &mode);    //parse input config file
    check_config(&config, mode);                //check that global variables were set correctly
    mk_run_nm(sample_name, &config);		    //construct sample name for shapemapper2
    
    if (mode == MULTI || mode == SINGLE) {
        print_MLT_SM2_script(sample_name, &config); //generate shapemapper2 run script
        
    } else if (mode == MULTIPLEX) {
        //parse_brcd_id_list();
        
    } else {
        printf("mk_run_script: error - unexpected run mode. aborting...\n");
        abort();
    }
    
    printf("%s\n\n", sample_name);    //print final sample name
    
    return 1;
}

/* check initialization: check initialization of config struct */
void check_cnfgMLT_init(tprobe_configuration * config)
{
    printf("\ninitial config values\n");
    printf("shapemapper2_path\t%2d\n",      config->shpmppr2_path[0]);
    
    printf("input_file_location\t%2d\n",    config->ipt_file_loc[0]);
    printf("output_file_location\t%2d\n",   config->out_file_loc[0]);
    
    printf("untreated_read1_prefix\t%2d\n", config->untreated_read1_prefix[0]);
    printf("untreated_read2_prefix\t%2d\n", config->untreated_read2_prefix[0]);
    
    printf("modified_read1_prefix\t%2d\n",  config->modified_read1_prefix[0]);
    printf("modified_read2_prefix\t%2d\n",  config->modified_read2_prefix[0]);
    
    printf("target_files_location\t%2d\n",  config->trg_files_loc[0]);
    printf("target_files_prefix\t%2d\n",    config->trg_files_prfx[0]);
    printf("min_target_length\t%2d\n",      config->min_target_length);
    printf("max_target_length\t%2d\n",      config->max_target_length);
    
    printf("input_name\t\t%2d\n",           config->input_name[0]);
    
    printf("cotranscriptional\t%2d\n",      config->cotranscriptional);
    printf("chemical_probe\t\t%2d\n",       config->chemical_probe[0]);
    
    printf("concatenated\t\t%2d\n",         config->concatenated);
    printf("runID\t\t\t%2d\n",              config->runID[0][0]);
    printf("smoothing\t\t%2d\n",            config->smoothing);
    
    printf("ligand_name\t\t%2d\n",          config->ligand_name[0]);
    printf("ligand_conc\t\t%2d\n",          config->ligand_conc[0]);
    
    printf("field\t\t\t%2d\n",              config->field[0][0]);
    printf("run_count\t\t%2d\n",            config->run_count);
    printf("field_count\t\t%2d\n\n",        config->field_count);
}

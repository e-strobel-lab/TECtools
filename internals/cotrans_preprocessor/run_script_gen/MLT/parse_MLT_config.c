//
//  parse_MLT_config.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/gen_utils.h"
#include "../../../utils/io_management.h"
#include "config_MLT_struct.h"

#include "parse_MLT_config.h"

/* parse_MLT_config: set config_MLT structure variables to config file settings */
int parse_MLT_config(FILE *ifp, configuration_MLT * config_MLT)
{
    int i = 0;
    int j = 0;
    
    char line[MAX_LINE] = {0};		//array for storing config file line
    char setting[MAX_LINE] = {0};	//array for storing config file setting name
    char val[MAX_LINE] = {0};		//array for storing config file setting value
    
    char * setting_start = {NULL};	//pointer to start of confige file setting in input line
        
    //in each iteration, get a line from the config file, check if it contains
    //a setting, and set the corresponding value in the config_MLT struct
    while (get_line(line, ifp)) {
        
        //reset setting and val variables
        setting[0] = '\0';
        val[0] = '\0';
        
        // > indicates data line in config file
        if (line[0] == '>') {
            
            //identify setting name in input line
            for (i = 1; isspace(line[i]) && line[i] && i < MAX_LINE; i++) {}; //skip to setting name
            if (!isspace(line[i]) && (i == 2 || i == 3)) { //setting name always starts at index 2 or 3
                setting_start = &line[i]; //set pointer to start of setting name
            } else {
                printf("parse_MLT_config: error unexpected setting line format. aborting... \n");
                abort();
            }
            
            //copy setting name to setting array
            for (i = 0; setting_start[i] != '=' && setting_start[i] && i < MAX_LINE; i++) {
                setting[i] = setting_start[i];
            }
            setting[i] = '\0';
            if (setting_start[i] != '=') { //check that loop ended on '=' character
                printf("parse_MLT_config: error unexpected setting name format. aborting... \n");
                abort();
            }
            
            //copy setting value to val array
            for (i++, j = 0; setting_start[i] && i < MAX_LINE; i++, j++) {
                val[j] = setting_start[i];
            }
            val[j] = '\0';
            if (setting_start[i]) { //check that loop ended on null character
                printf("parse_MLT_config: error unexpected setting value format. aborting... \n");
                abort();
            }
            
            //check for setting value errors
            //required settings have a default value of "FILL_IN"
            //optional settings have a default value of "NULL"
            //and no setting should have an empty value
            if (!strcmp(val, "FILL_IN")) { //check if required setting value has been omitted
                printf("parse_MLT_config: error - missing required value for %s setting. aborting...\n", setting);
                abort();
            } else if (val[0] == '\0') {   //check of setting value is empty in config file
                printf("parse_MLT_config: error - value for %s setting is empty. optional settings without a value must be set to NULL. aborting...\n", setting);
                abort();
            } else if (!strcmp(val, "NULL")) {
                //nothing to do here, just forcing NULL values away from the if-else statement below
            } else {
                
                //set path to shapemapper2 executable
                if (!strcmp(setting, "shapemapper2_path")) {
                    strcpy(config_MLT->shpmppr2_path, val);
                    
                //set path to input files
                } else if (!strcmp(setting, "input_file_location")) {
                    strcpy(config_MLT->ipt_file_loc, val);
                    chk_filepath_frmt(config_MLT->ipt_file_loc);
                    
                //set path to output file location
                } else if (!strcmp(setting, "output_file_location")) {
                    strcpy(config_MLT->out_file_loc, val);
                    chk_filepath_frmt(config_MLT->out_file_loc);
                    
                //set prefix for untreated read 1 files
                } else if (!strcmp(setting, "untreated_read1_prefix")) {
                    strcpy(config_MLT->untreated_read1_prefix, val);
                    
                //set prefix for untreated read 2 files
                } else if (!strcmp(setting, "untreated_read2_prefix")) {
                    strcpy(config_MLT->untreated_read2_prefix, val);
                    
                //set prefix for modified read 1 files
                } else if (!strcmp(setting, "modified_read1_prefix")) {
                    strcpy(config_MLT->modified_read1_prefix, val);
                    
                //set prefix for modified read 2 files
                } else if (!strcmp(setting, "modified_read2_prefix")) {
                    strcpy(config_MLT->modified_read2_prefix, val);
                    
                //set path to target files
                } else if (!strcmp(setting, "target_files_location")) {
                    strcpy(config_MLT->trg_files_loc, val);
                    chk_filepath_frmt(config_MLT->trg_files_loc);
                    
                //set target file prefix
                } else if (!strcmp(setting, "target_files_prefix")) {
                    strcpy(config_MLT->trg_files_prfx, val);
                    
                //set minimum observed target length
                } else if (!strcmp(setting, "min_target_length")) {
                    config_MLT->min_target_length = atoi(val);
                    
                //set maximum target length
                } else if (!strcmp(setting, "max_target_length")) {
                    config_MLT->max_target_length = atoi(val);
                    
                //set input name
                } else if (!strcmp(setting, "name")) {
                    strcpy(config_MLT->input_name, val);
                    
                //set cotranscriptional flag
                } else if (!strcmp(setting, "cotranscriptional")) {
                    set_TF_value(val, setting, &(config_MLT->cotranscriptional));

                //set chemical probe identity
                } else if (!strcmp(setting, "chemical_probe")) {
                    strcpy(config_MLT->chemical_probe, val);
                    
                //set concatenated flag
                } else if (!strcmp(setting, "concatenated")) {
                    set_TF_value(val, setting, &(config_MLT->concatenated));
                    
                //set runID values
                } else if (!strcmp(setting, "runID")) {
                    if (strcmp(val, "NULL")) {                  //if value is not NULL
                        if (config_MLT->run_count < MAX_RUNS) { //and if MAX_RUNS has not yet been reached
                            if (strstr(val, "_") == NULL) {     //and if the runID does not contain an underscore
                                strcpy(config_MLT->runID[config_MLT->run_count++], val); //store runID value
                            } else { //underscore detected, abort
                                printf("parse_MLT_config: error - runID cannot contain an underscore character. aborting...\n");
                                abort();
                            }
                        } else { //too many runIDs, abort
                            printf("parse_MLT_config: error - runID count exceeds %d. aborting...\n", MAX_RUNS);
                            abort();
                        }
                    }
                    
                //set smoothing flag
                } else if (!strcmp(setting, "smoothing")) {
                    set_TF_value(val, setting, &(config_MLT->smoothing));
                    
                //set ligand name (if ligand was used)
                } else if (!strcmp(setting, "ligand_name")) {
                    strcpy(config_MLT->ligand_name, val);
                    
                //set ligand concentration (if ligand was used)
                } else if (!strcmp(setting, "ligand_conc")) {
                    strcpy(config_MLT->ligand_conc, val);
                    
                //set optional custom field values
                } else if (!strcmp(setting, "field")) {
                    if (strcmp(val, "NULL")) {                      //if value is not NULL
                        if (config_MLT->field_count < MAX_FIELDS) { //and if MAX_FIELDS has not yet been reached
                            if (strstr(val, "_") == NULL) {         //and if the field does not contain a uscore
                                strcpy(config_MLT->field[config_MLT->field_count++], val);
                            } else { //underscore was detected, abort
                                printf("parse_MLT_config: error - custom field values cannot contain an underscore character. aborting...\n");
                                abort();
                            }
                        } else { //too many fields, abort
                            printf("parse_MLT_config: error - field count exceeds 8. aborting...\n");
                            abort();
                        }
                    }
                    
                //throw error for unexpected setting name
                } else {
                    printf("parse_MLT_config: error - unrecognized setting (%s). aborting...\n", setting);
                    abort();
                }
            }
        } else if (line[0] != '#') {
            printf("parse_MLT_config: error - unexpected line initator (%c). lines must start with # (comment) or > (data). aborting...\n", line[0]);
            abort();
        }
    }

    return 1;
}

/* set_TF_value: set true or false value as 1 and 0, respectively */
void set_TF_value(char * tf, char * setting_name, int * config_val)
{
    if (!strcmp(tf, "FALSE")) {
        *config_val = 0;
    } else if (!strcmp(tf, "TRUE")) {
        *config_val = 1;
    } else {
        printf("parse_MLT_config: error - unexpected value for %s setting. set to FALSE, or TRUE. aborting...\n", setting_name);
        abort();
    }
}









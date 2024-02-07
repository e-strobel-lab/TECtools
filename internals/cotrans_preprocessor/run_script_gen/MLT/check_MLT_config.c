//
//  check_MLT_config.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "config_MLT_struct.h"

#include "check_MLT_config.h"

//check_MLT_config: check that config_MLT variables have been set properly
int check_MLT_config(configuration_MLT * config_MLT)
{
    int i = 0;
    
    FILE * out_fp = NULL;
    
    //generate output parsed config file
    if ((out_fp = fopen("parsed_config.txt", "w")) == NULL) {
        printf("check_MLT_config: ERROR - could not generate parsed_config.txt file. Aborting program...\n");
        abort();
    }
    
    //check shpmppr2_path variable
    if (!config_MLT->shpmppr2_path[0]) {
        printf("check_MLT_config: error - shpmppr2_path is not set. aborting\n");
        abort();
    } else {
        printf("\nshapemapper2_path\t%s\n\n", config_MLT->shpmppr2_path);
        fprintf(out_fp, "\nshapemapper2_path\t%s\n\n", config_MLT->shpmppr2_path);
    }
    
    //check ipt_file_loc variable
    if (!config_MLT->ipt_file_loc[0]) {
        printf("check_MLT_config: error - ipt_file_loc is not set. aborting\n");
        abort();
    } else {
        printf("input_file_location\t%s\n", config_MLT->ipt_file_loc);
        fprintf(out_fp, "input_file_location\t%s\n", config_MLT->ipt_file_loc);
    }
    
    //check out_file_loc variable
    if (!config_MLT->out_file_loc[0]) {
        printf("check_MLT_config: error - out_file_loc is not set. aborting\n");
        abort();
    } else {
        printf("output_file_location\t%s\n\n", config_MLT->out_file_loc);
        fprintf(out_fp, "output_file_location\t%s\n\n", config_MLT->out_file_loc);
    }
    
    //check untreated_read1_prefix
    if (!config_MLT->untreated_read1_prefix[0]) {
        printf("check_MLT_config: error - untreated_read1_prefix is not set. aborting\n");
        abort();
    } else {
        printf("untreated_read1_prefix\t%s\n", config_MLT->untreated_read1_prefix);
        fprintf(out_fp, "untreated_read1_prefix\t%s\n", config_MLT->untreated_read1_prefix);
    }
    
    //check untreated_read2_prefix
    if (!config_MLT->untreated_read2_prefix[0]) {
        printf("check_MLT_config: error - untreated_read2_prefix is not set. aborting\n");
        abort();
    } else {
        printf("untreated_read2_prefix\t%s\n", config_MLT->untreated_read2_prefix);
        fprintf(out_fp, "untreated_read2_prefix\t%s\n", config_MLT->untreated_read2_prefix);
    }
    
    //check modified_read1_prefix
    if (!config_MLT->modified_read1_prefix[0]) {
        printf("check_MLT_config: error - modified_read1_prefix is not set. aborting\n");
        abort();
    } else {
        printf(" modified_read1_prefix\t%s\n", config_MLT->modified_read1_prefix);
        fprintf(out_fp, " modified_read1_prefix\t%s\n", config_MLT->modified_read1_prefix);
    }
    
    //check modified_read2_prefix
    if (!config_MLT->modified_read2_prefix[0]) {
        printf("check_MLT_config: error - modified_read2_prefix is not set. aborting\n");
        abort();
    } else {
        printf(" modified_read2_prefix\t%s\n\n", config_MLT->modified_read2_prefix);
        fprintf(out_fp, " modified_read2_prefix\t%s\n\n", config_MLT->modified_read2_prefix);
    }
    
    //check trg_files_loc
    if (!config_MLT->trg_files_loc[0]) {
        printf("check_MLT_config: error - trg_files_loc is not set. aborting\n");
        abort();
    } else {
        printf("target_files_location\t%s\n", config_MLT->trg_files_loc);
        fprintf(out_fp, "target_files_location\t%s\n", config_MLT->trg_files_loc);
    }
    
    //check trg_files_prfx
    if (!config_MLT->trg_files_prfx[0]) {
        printf("check_MLT_config: error - trg_files_prfx is not set. aborting\n");
        abort();
    } else {
        printf("target_files_prefix\t%s\n", config_MLT->trg_files_prfx);
        fprintf(out_fp, "target_files_prefix\t%s\n", config_MLT->trg_files_prfx);
    }
    
    //check min_target_length
    if (config_MLT->min_target_length == -1) {
        printf("check_MLT_config: error - min_target_length is not set. aborting\n");
        abort();
    } else {
        printf("min_target_length\t%03d\n", config_MLT->min_target_length);
        fprintf(out_fp, "min_target_length\t%03d\n", config_MLT->min_target_length);
    }
    
    //check max_target_length
    if (config_MLT->max_target_length == -1) {
        printf("check_MLT_config: error - max_target_length is not set. aborting\n");
        abort();
    } else {
        printf("max_target_length\t%03d\n\n", config_MLT->max_target_length);
        fprintf(out_fp, "max_target_length\t%03d\n\n", config_MLT->max_target_length);
    }
    
    //check input_name
    if (!config_MLT->input_name[0]) {
        printf("check_MLT_config: error - input_name is not set. aborting\n");
        abort();
    } else {
        printf("input_name\t\t%s\n", config_MLT->input_name);
        fprintf(out_fp, "input_name\t\t%s\n", config_MLT->input_name);
    }
    
    //check cotranscriptional
    if (config_MLT->cotranscriptional == -1) {
        printf("check_MLT_config: error - cotranscriptional is not set. aborting\n");
        abort();
    } else {
        printf("cotranscriptional\t");
        fprintf(out_fp, "cotranscriptional\t");
        switch (config_MLT->cotranscriptional) {
            case  0:
                printf("FALSE\n");
                fprintf(out_fp, "FALSE\n");
                break;
            case  1:
                printf("TRUE\n");
                fprintf(out_fp, "TRUE\n");
                break;
            default: break;
        }
    }
    
    //check chemical_probe
    if (!config_MLT->chemical_probe[0]) {
        printf("check_MLT_config: error - chemical_probe is not set. aborting\n");
        abort();
    } else {
        printf("chemical_probe\t\t%s\n\n", config_MLT->chemical_probe);
        fprintf(out_fp, "chemical_probe\t\t%s\n\n", config_MLT->chemical_probe);
    }
    
    //check concatenated
    if (config_MLT->concatenated == -1) {
        printf("check_MLT_config: error - concatenated is not set. aborting\n");
        abort();
    } else {
        if (config_MLT->concatenated == 1 && config_MLT->run_count < 2) {
            printf("check_MLT_config: error - concatenated=TRUE but <2 runIDs were provided. aborting...\n");
            abort();
        } else if (!config_MLT->concatenated && config_MLT->run_count > 1) {
            printf("check_MLT_config: error - concatenated=FALSE but >1 runID was provided. aborting...\n");
            abort();
        } else {
            printf("concatenated\t\t");
            fprintf(out_fp, "concatenated\t\t");
            switch (config_MLT->concatenated) {
                case -1:
                    printf("NULL\n");
                    fprintf(out_fp, "NULL\n");
                    break;
                case  0:
                    printf("FALSE\n");
                    fprintf(out_fp, "FALSE\n");
                    break;
                case  1:
                    printf("TRUE\n");
                    fprintf(out_fp, "TRUE\n");
                    break;
                default: break;
            }
        }
    }
    
    //check runID
    if (!config_MLT->runID[0][0]) {
        printf("check_MLT_config: error - runID is not set. aborting...\n");
        abort();
    } else {
        printf("run_count\t\t%d\n", config_MLT->run_count);
        fprintf(out_fp, "run_count\t\t%d\n", config_MLT->run_count);
        for (i = 0; i < config_MLT->run_count; i++) {
            printf("runID%d\t\t\t%s\n", i, config_MLT->runID[i]);
            fprintf(out_fp, "runID%d\t\t\t%s\n", i, config_MLT->runID[i]);
        }
    }
    
    //check ligand_name and ligand_conc settings
    if (config_MLT->ligand_name[0] && !config_MLT->ligand_conc[0]) {
        printf("check_MLT_config: error - ligand name was supplied without ligand concentrtation. aborting...\n");
        abort();
    } else if (!config_MLT->ligand_name[0] && config_MLT->ligand_conc[0]) {
        printf("check_MLT_config: error - ligand concentration was supplied without ligand name. aborting...\n");
        abort();
	//check format of ligand concentration value
    } else if (!isdigit(config_MLT->ligand_conc[0]) || //first three conc chars must be digits
               !isdigit(config_MLT->ligand_conc[1]) ||
               !isdigit(config_MLT->ligand_conc[2]) ||
               
               //molarity prefix must be m, u, n, or p
               (config_MLT->ligand_conc[3] != 'm'&&
                config_MLT->ligand_conc[3] != 'u'&&
                config_MLT->ligand_conc[3] != 'n'&&
                config_MLT->ligand_conc[3] != 'p')  ||
            	
               //conc unit must be M
               config_MLT->ligand_conc[4] != 'M') {
        
        printf("check_MLT_config: error - unexpected ligand concentration format. aborting...\n");
        abort();
    } else if (config_MLT->ligand_name[0] && config_MLT->ligand_conc[0]) {
        printf("ligand_name\t\t%s\n", config_MLT->ligand_name);
        fprintf(out_fp, "ligand_name\t\t%s\n", config_MLT->ligand_name);
        
        printf("ligand_conc\t\t%s\n\n", config_MLT->ligand_conc);
        fprintf(out_fp, "ligand_conc\t\t%s\n\n", config_MLT->ligand_conc);
    }
    
    //check smoothing
    if (config_MLT->smoothing == -1) {
        printf("check_MLT_config: error - smoothing is not set. aborting...\n");
        abort();
    } else {
        printf("smoothing\t\t");
        fprintf(out_fp, "smoothing\t\t");
        switch (config_MLT->smoothing) {
            case  0:
                printf("FALSE\n\n");
                fprintf(out_fp, "FALSE\n\n");
                break;
            case  1:
                printf("TRUE\n\n");
                fprintf(out_fp, "TRUE\n\n");
                break;
            default: break;
        }
    }
    
    //print field settings
    for (i = 0; i < config_MLT->field_count; i++) {
        printf("field%d\t\t\t%s\n", i, config_MLT->field[i]);
        fprintf(out_fp, "field%d\t\t\t%s\n", i, config_MLT->field[i]);
    }
    
    if (i) { //printed field settings
        printf("\n");
        fprintf(out_fp, "\n");
    }
    
    //close output file
    if ((fclose(out_fp)) == EOF) {
        printf("check_MLT_config: ERROR - error occurred when closing parsed_config.txt file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

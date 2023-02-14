//
//  print_MLT_SM2_script.c
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
#include "config_MLT_struct.h"

#include "print_MLT_SM2_script.h"

/*print_MLT_SM2_script: print shell script that runs shapemapper2 for all 3' end lengths */
int print_MLT_SM2_script(char * sample_name, configuration_MLT * config_MLT)
{
    int i = 0;
    
    char crrnt_UNT_READ1[MAX_LINE*2] = {0};	//array for current untreated read 1 filename
    char crrnt_UNT_READ2[MAX_LINE*2] = {0};	//array for current untreated read 2 filename
    char crrnt_MOD_READ1[MAX_LINE*2] = {0};	//array for current modified read 1 filename
    char crrnt_MOD_READ2[MAX_LINE*2] = {0};	//array for current modified read 2 filename
    
    char * suffix = NULL;						 //pointer to fastq file suffix
    char with_smoothing[MAX_LINE] = "_SM.fq.gz";	 //suffix for fastq files with smoothing
    char without_smoothing[MAX_LINE] = ".fq.gz"; //suffix for fastq files without smoothing
    
    //set suffix pointer based on smoothing option in config file
    if (config_MLT->smoothing) {
        suffix = &with_smoothing[0];
    } else {
        suffix = &without_smoothing[0];
    }
    
    FILE * out_fp = NULL; //output file pointer
    
    char script_name[MAX_LINE] = {0};				//array for run script name
    sprintf(script_name, "%s_run.sh", sample_name); //run script is named using sample name
    
    //generate output shell script file
    if ((out_fp = fopen(script_name, "w")) == NULL) {
        printf("print_MLT_SM2_script: ERROR - could not generate config file. Aborting program...\n");
        abort();
    }
    
    //write output shell script file
    fprintf(out_fp, "#!/bin/bash\n"); //print she-bang
    
    //print shapemapper2 run command for each 3' end from min to max target length
    for (i = config_MLT->min_target_length; i <= config_MLT->max_target_length; i++) {
        
        //populate file name strings for current loop iteration
        sprintf(crrnt_UNT_READ1, "%s_%03d_R1%s", config_MLT->untreated_read1_prefix, i, suffix);
        sprintf(crrnt_UNT_READ2, "%s_%03d_R2%s", config_MLT->untreated_read2_prefix, i, suffix);
        sprintf(crrnt_MOD_READ1, "%s_%03d_R1%s", config_MLT->modified_read1_prefix,  i, suffix);
        sprintf(crrnt_MOD_READ2, "%s_%03d_R2%s", config_MLT->modified_read2_prefix,  i, suffix);
        
        //print directory commands to output shell script
        fprintf(out_fp, "mkdir %s%03d_analysis\n", config_MLT->out_file_loc, i);
        fprintf(out_fp, "cd %s%03d_analysis\n",    config_MLT->out_file_loc, i);
        
        /***** print shapemapper 2 run command *****/
        
        //shapemapper 2 path
        fprintf(out_fp, "%s", config_MLT->shpmppr2_path);
        
        //sample name
        fprintf(out_fp, " --name %s_%03d", sample_name, i);
        
        //target file name
        fprintf(out_fp, " --target %s%s_%03dnt_target.fa",
                config_MLT->trg_files_loc, config_MLT->trg_files_prfx, i);
        
        //output name
        fprintf(out_fp, " --out %s_%03d_out", sample_name, i);
        
        //minimum depth option
        fprintf(out_fp, " --min-depth 500");
        
        //minimum quality to trim option
        fprintf(out_fp, " --min-qual-to-trim 10");
        
        //minimum quality to count option
        fprintf(out_fp, " --min-qual-to-count 25");
        
        //modified read inputs
        fprintf(out_fp, " --modified");
        fprintf(out_fp, " --R1 %s%s", config_MLT->ipt_file_loc, crrnt_MOD_READ1);
        fprintf(out_fp, " --R2 %s%s", config_MLT->ipt_file_loc, crrnt_MOD_READ2);
        
        //untreated read inputs
        fprintf(out_fp, " --untreated");
        fprintf(out_fp, " --R1 %s%s", config_MLT->ipt_file_loc, crrnt_UNT_READ1);
        fprintf(out_fp, " --R2 %s%s", config_MLT->ipt_file_loc, crrnt_UNT_READ2);
        
        //send screen prints to outfile
        fprintf(out_fp, " > outfile_%03d.txt", i);
        
        //run command in background, end command
        fprintf(out_fp, " &\n");
    }
    
    if ((fclose(out_fp)) == EOF) {
        printf("print_MLT_SM2_script: ERROR - error occurred when closing shell script file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

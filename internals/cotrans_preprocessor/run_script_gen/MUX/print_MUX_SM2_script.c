//
//  print_MUX_SM2_script.c
//  
//
//  Created by Eric Strobel on 9/24/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../UNV/config_struct.h"

#include "print_MUX_SM2_script.h"


/*print_MUX_SM2_script: print shell script that runs shapemapper2 for all barcoded targets */
int print_MUX_SM2_script(char * sample_name, int brcd_cnt, char ** brcd_id, tprobe_configuration * config)
{
    int i = 0;   //general purpose index
    int l = 0;   //layer index
    int ret = 0; //snprintf return value
    
    int cmnds_wrtn = 0; //number of SM2 run commands written
    
    char crrnt_UNT_READ1[MAX_LINE*2] = {0}; //array for current untreated read 1 filename
    char crrnt_UNT_READ2[MAX_LINE*2] = {0}; //array for current untreated read 2 filename
    char crrnt_MOD_READ1[MAX_LINE*2] = {0}; //array for current modified read 1 filename
    char crrnt_MOD_READ2[MAX_LINE*2] = {0}; //array for current modified read 2 filename
    
    char run_script_nm[MAX_LINE+1] = {0}; //array to store run script name
    char suffix[7] = ".fq.gz";            //suffix for fastq files
    
    int parallel = 100;                                                   //num of SM2 runs to run in parallel
    int layers = (brcd_cnt / parallel) + ((brcd_cnt % parallel) ? 1 : 0); //num of commands in each script
        
    FILE  * coord_fp = NULL; //output file pointer for coordination shell script
    FILE ** out_fp = NULL;   //output file pointer for run SM2 shell scripts
    
    //open coordination shell script file
    if ((coord_fp = fopen("run_SM2.sh", "w")) == NULL) {
        printf("print_MUX_SM2_script: error - failed to generate run_SM2.sh file. aborting...\n");
        abort();
    }
    fprintf(coord_fp, "#!/bin/bash\n"); //print she-bang

    //allocate file pointers for run scripts
    if ((out_fp = calloc(parallel, sizeof(*out_fp))) == NULL) {
        printf("print_MUX_SM2_script: error - memory allocation for output file pointers failed. aborting...\n");
        abort();
    }
    
    cmnds_wrtn = 0; //initialize commands written counter
    
    //open run script files and write coordination shell script
    for (i = 0; i < parallel; i++) {
        
        //generate run script file name
        ret = snprintf(run_script_nm, MAX_LINE, "script_%03d.sh", i+1);
        if (ret >= MAX_LINE || ret < 0) {
            printf("print_MUX_SM2_script: error - error when generating run script name. aborting...\n");
            abort();
        }
         
        //open run script file
        if ((out_fp[i] = fopen(run_script_nm, "w")) == NULL) {
            printf("print_MUX_SM2_script: error - failed to generate SM2 run script file. aborting...\n");
            abort();
        }
        
        fprintf(out_fp[i], "#!/bin/bash\n");           //print she-bang
        fprintf(coord_fp, "sh %s &\n", run_script_nm); //write run command to coordination shell script
    }
    
    //write run scripts
    for (l = 0; l < layers && cmnds_wrtn < brcd_cnt; l++) {
        
        for (i = 0; i < parallel && cmnds_wrtn < brcd_cnt; i++) {
            
            //populate file name strings for current loop iteration
            //TODO: add safeguard to parse_brcd_id_list that ensures barcode ids are 5 chars long?
            sprintf(crrnt_UNT_READ1, "%s_%s_R1%s", config->untreated_read1_prefix, brcd_id[cmnds_wrtn], suffix);
            sprintf(crrnt_UNT_READ2, "%s_%s_R2%s", config->untreated_read2_prefix, brcd_id[cmnds_wrtn], suffix);
            sprintf(crrnt_MOD_READ1, "%s_%s_R1%s", config->modified_read1_prefix,  brcd_id[cmnds_wrtn], suffix);
            sprintf(crrnt_MOD_READ2, "%s_%s_R2%s", config->modified_read2_prefix,  brcd_id[cmnds_wrtn], suffix);
            
            //print directory commands to output shell script
            fprintf(out_fp[i], "mkdir %s%s_analysis\n", config->out_file_loc, brcd_id[cmnds_wrtn]);
            fprintf(out_fp[i], "cd %s%s_analysis\n",    config->out_file_loc, brcd_id[cmnds_wrtn]);
            
            //print shapemapper 2 run command
            
            //shapemapper 2 path
            fprintf(out_fp[i], "%s", config->shpmppr2_path);
            
            //sample name
            fprintf(out_fp[i], " --name %s_%s", sample_name, brcd_id[cmnds_wrtn]);
            
            //target file name
            fprintf(out_fp[i], " --target %s%s_%s.fa",
                    config->trg_files_loc, config->trg_files_prfx, brcd_id[cmnds_wrtn]);
            
            //output name
            fprintf(out_fp[i], " --out %s_%s_out", sample_name, brcd_id[cmnds_wrtn]);
            
            //minimum depth option
            fprintf(out_fp[i], " --min-depth 500");
            
            //minimum quality to trim option
            fprintf(out_fp[i], " --min-qual-to-trim 10");
            
            //minimum quality to count option
            fprintf(out_fp[i], " --min-qual-to-count 25");
            
            //modified read inputs
            fprintf(out_fp[i], " --modified");
            fprintf(out_fp[i], " --R1 %s%s", config->ipt_file_loc, crrnt_MOD_READ1);
            fprintf(out_fp[i], " --R2 %s%s", config->ipt_file_loc, crrnt_MOD_READ2);
            
            //untreated read inputs
            fprintf(out_fp[i], " --untreated");
            fprintf(out_fp[i], " --R1 %s%s", config->ipt_file_loc, crrnt_UNT_READ1);
            fprintf(out_fp[i], " --R2 %s%s", config->ipt_file_loc, crrnt_UNT_READ2);
            
            //send screen prints to outfile
            fprintf(out_fp[i], " > outfile_%s.txt\n", brcd_id[cmnds_wrtn]);
            
            cmnds_wrtn++;
        }
    }
    
    //close run scripts
    for (i = 0; i < parallel; i++) {
        if (fclose(out_fp[i]) == EOF) {
            printf("print_MUX_SM2_script: error - failed to close SM2 run script file. aborting...\n");
            abort();
        }
    }

    //close coordination shell script
    if ((fclose(coord_fp)) == EOF) {
        printf("print_MLT_SM2_script: ERROR - error occurred when closing run_SM2.sh shell script file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

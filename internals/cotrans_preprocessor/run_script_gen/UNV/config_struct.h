//
//  config_struct.h
//  
//
//  Created by Eric Strobel on 3/18/22.
//

#ifndef config_struct_h
#define config_struct_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#define MAX_RUNS 8		//maximum sequencing run IDs
#define MAX_FIELDS 8	//maximum custom fields

typedef struct tprobe_configuration {
    char * shpmppr2_path;          //path to shapemapper 2 executable
    
    char * ipt_file_loc;           //path to directory that contains input files
    char * out_file_loc;           //path to directory that will be used for output files
    
    char * untreated_read1_prefix; //prefix for untreated read 1 fastq files
    char * untreated_read2_prefix; //prefix for untreated read 2 fastq files
    char * modified_read1_prefix;  //prefix for modified read 1 fastq files
    char * modified_read2_prefix;  //prefix for modified read 2 fastq files
    
    char * trg_files_loc;          //path to directory that contains target files
    char * trg_files_prfx;         //prefix for target files
    int min_target_length;         //minimum target length
    int max_target_length;         //maximum target length
    
    char * input_name;             //user-specified name for sample
    int  cotranscriptional;		   //flag indicating cotranscriptional or equilibrium structure probing
    char * chemical_probe;         //array for chemical probe identifier
    
    int  concatenated;             //flag that indicates if fastq files were concatenated from multiple seq runs
    char * runID[MAX_RUNS];        //identifiers for the sequencing run(s) that generated the fastq files
    int  smoothing;		           //flag that indicates if the analysis will be performed using smoothed fastqs
    char * ligand_name;	           //array for storing ligand name
    char * ligand_conc;            //array for storing ligand concentration
    char * field[MAX_FIELDS];      //arrays for optional user-specified sample name information
    
    int run_count;	 //number of sequencing runIDs, can be >1 if reads were concatenated from multiple runs
    int field_count; //number of user-specified fields
    
    
} tprobe_configuration;

/* set_cfg_string: allocate memory for configuration_MLT structure string */
void set_cfg_string(char ** cfg_str, char * val, int buffer);

/* set_TF_value: set true or false value as 1 and 0, respectively */
void set_TF_value(char * tf, char * setting_name, int * config_val);

#endif /* config_struct_h */

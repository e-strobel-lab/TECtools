//
//  config_MLT_struct.h
//  
//
//  Created by Eric Strobel on 3/18/22.
//

#ifndef config_MLT_struct_h
#define config_MLT_struct_h

#include <stdio.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#define MAX_RUNS 8		//maximum sequencing run IDs
#define MAX_FIELDS 8	//maximum custom fields

typedef struct configuration_MLT {
    char shpmppr2_path[MAX_LINE];			//path to shapemapper 2 executable
    
    char ipt_file_loc[MAX_LINE];			//path to directory that contains input files
    char out_file_loc[MAX_LINE];			//path to directory that will be used for output files
    
    char untreated_read1_prefix[MAX_LINE];	//prefix for untreated read 1 fastq files
    char untreated_read2_prefix[MAX_LINE];	//prefix for untreated read 2 fastq files
    char modified_read1_prefix[MAX_LINE];	//prefix for modified read 1 fastq files
    char modified_read2_prefix[MAX_LINE];	//prefix for modified read 2 fastq files
    
    char trg_files_loc[MAX_LINE];			//path to directory that contains target files
    char trg_files_prfx[MAX_LINE];			//prefix for target files
    int  min_target_length;					//minimum target length
    int  max_target_length;					//maximum target length
    
    char input_name[MAX_LINE];	   //user-specified name for sample
    int  cotranscriptional;		   //flag indicating cotranscriptional or equilibrium structure probing
    char chemical_probe[MAX_LINE]; //array for chemical probe identifier
    
    int  concatenated;	//flag that indicates if fastq files were concatenated from multiple seq runs
    char runID[MAX_RUNS][MAX_LINE];	//identifiers for the sequencing run(s) that generated the fastq files
    int  smoothing;		//flag that indicates if the analysis will be performed using smoothed fastqs
    char ligand_name[MAX_LINE];	//array for storing ligand name
    char ligand_conc[MAX_LINE];	//array for storing ligand concentration
    char field[MAX_FIELDS][MAX_LINE]; //arrays for optional user-specified sample name information
    
    int run_count;	 //number of sequencing runIDs, can be >1 if reads were concatenated from multiple runs
    int field_count; //number of user-specified fields
    
    
} configuration_MLT;

#endif /* config_MLT_struct_h */

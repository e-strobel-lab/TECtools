//
//  process_TECprobeVL_profiles.c
//  
//
//  Created by Eric Strobel on 1/19/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include <dirent.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"
#include "../../utils/io_management.h"

#include "../../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "../global/store_SM2_profile.h"
#include "../global/initialize_empty_profile.h"
#include "../global/calculate_normalization_factor.h"

#include "./parse_VL_sample_name.h"
#include "./read_VL_analysis_directories.h"
#include "./VL_input_validation.h"
#include "./make_VL_output_directories.h"
#include "./normalize_VL_reactivities.h"
#include "./merge_VL_profiles.h"
#include "./generate_VL_sample_name.h"
#include "./print_merged_VL_profiles.h"
#include "./print_legacy_compiled_table.h"

void print_processing_record(sample_names * sn, output_files * outfiles, SM2_analysis_directory * an_dir, SM2_analysis_directory * mrg);

int main(int argc, char *argv[])
{
    extern int debug; //debug mode flag
    
    char * prnt_dir_nm[MAX_RUNS] = {NULL}; //list of parent directories
    
    int dir_count = 0; //flag that input directory was supplied
    
    int min_depth = 1000;  //minimum effective read depth
                           //NOTE: SM2 docs say default 5000, but 1000 is used in code
    double max_bkg = 0.05; //maximum untreated mutation rate
    
    int min_depth_provided = 0; //flag that minimum depth option was provided
    int max_bkg_provided = 0;   //flag that maximum background option was provided
    int verify_norm = 0;        //flag to verify normalization against SM2 calculations
    
    sample_names sn = {{{0}}};  //structure for sample name storage/merged name construction
    output_files outfiles;      //structure for storing output file pointers and names
    
    char dflt_out_dir_nm[20] = {"dataset_norm_out"};
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"input",          required_argument, 0, 'i'}, //shapemapper analysis directory
            {"out_dir_name",   required_argument, 0, 'o'}, //output directory name
            {"sample-name",    required_argument, 0, 'n'}, //user-supplied sample name
            {"debug",          no_argument,       0, 'd'}, //debug flag
            {"min-depth",      required_argument, 0, 'e'}, //min effective depth
            {"max-background", required_argument, 0, 'b'}, //max background
            {"verify_norm",    no_argument,       0, 'y'}, //verify normalization
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "i:o:n:de:b:y", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
            case 'i': //input directory
                if (dir_count < MAX_RUNS) { //if the number of input directories has not exceed the max
                    
                    //allocate memory for parent directory name storage
                    if ((prnt_dir_nm[dir_count] = calloc(strlen(argv[optind-1])+1, sizeof(*(prnt_dir_nm[dir_count])))) == NULL) {
                        printf("process_TECprobeVL_profiles: error - memory allocation for directory name failed. aborting...\n");
                        abort();
                    }
                    strcpy(prnt_dir_nm[dir_count], argv[optind-1]); //store input directory name
                    dir_count++; //increment directory count
                    
                } else { //max input directories exceeded
                    printf("process_TECprobeVL_profiles: error - too many input directories were provided. The maximum number allowed is %d. aborting...\n", MAX_RUNS);
                    abort();
                }
                
                break;
                
            case 'o': //output directory name
                
                //store output directory name
                if (snprintf(outfiles.out_dir, MAX_NAME, "%s_%s", argv[optind-1], dflt_out_dir_nm) >= MAX_NAME) {
                    printf("process_TECprobeVL_profiles: error - output directory name exceeds buffer. aborting...\n");
                    abort();
                }
                
                break;
                
            case 'n': //user-supplied output file name
                
                //store user-supplied output file name
                if (snprintf(sn.usr, MAX_NAME, "%s", argv[optind-1]) >= MAX_NAME) {
                    printf("process_TECprobeVL_profiles: error - sample name exceeds buffer. aborting...\n");
                    abort();
                }
                break;
            
            case 'd': //turn on debug mode
                debug = 1;
                break;
            
            case 'e': //set minimum effective depth cutoff
                if (!min_depth_provided) {
                    check_int_str(argv[optind-1], ABORT_FAILURE);
                    min_depth = atoi(argv[optind-1]);
                } else {
                    printf("process_TECprobeVL_profiles: error - min-depth option was provided more than once. aborting...\n");
                    abort();
                }
                
                
                break;
                
            case 'b': //set maximum background mutation rate cutoff
                if (!max_bkg_provided) {
                    check_float_str(argv[optind-1], ABORT_FAILURE);
                    max_bkg = strtof(argv[optind-1], NULL);
                } else {
                    printf("process_TECprobeVL_profiles: error - max-background option was provided more than once. aborting...\n");
                    abort();
                }
                
                break;
                
            case 'y': //turn on normalization verification
                verify_norm = 1;
                break;
            
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    
    if (optind < argc) {
        printf("\nnon-option ARGV-elements:\n");
        while (optind < argc)
            printf("\n%s \n", argv[optind++]);
        putchar('\n');
        printf("Aborting program.\n\n");
        abort();
    }
    /*********** end of option parsing ***********/
    
    //check that input directory was provided
    if (dir_count < 1) {
        printf("process_TECprobeVL_profiles: error - too few input directories (%d) were provided. aborting...\n" ,dir_count);
        abort();
    }
    
    //if no output directory name as provided, set the output directory name
    if (!outfiles.out_dir[0]) {
        strcpy(outfiles.out_dir, dflt_out_dir_nm);
    }
    
    //set the output sample name that will be used
    if (sn.usr[0]) {            //if a user-supplied output sample name was provided
        sn.sn2use = &sn.usr[0]; //use the user-supplied output sample name
    } else {
        sn.sn2use = &sn.mrg[0]; //otherwise, use the auto-generated sample name
    }
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
        
    SM2_analysis_directory * an_dir = NULL; //pointer for SM2 analysis directory memory allocation
    
    //allocate memory for storing SM2 analysis directory information
    if ((an_dir = calloc(dir_count, sizeof(*an_dir))) == NULL) {
        printf("process_TECprobeVL_profiles: error - memory allocation for directory/file pointers failed. aborting...\n");
        abort();
    }
    
    //initialize SM2 analysis directory structs
    for (i = 0; i < dir_count; i++) {                  //for each input directory
        strcpy(an_dir[i].prnt_dir_nm, prnt_dir_nm[i]); //store the parent directory name
        an_dir[i].min_tl = MIN_INIT;                   //initialize min_tl to MIN_INIT (512 (MAX_ROW) + 1)
        an_dir[i].max_tl = MAX_INIT;                   //initialize max_tl to MAX_INIT (0)
    }
        
    //read each SM2 analysis directory and confirm that the number of
    //files and transcript length range is the same for all samples
    for (i = 0; i < dir_count; i++) {
        
        //open all relevant directories and files in the parent directory
        //then validate that transcript length profiles are contiguous
        an_dir[i].prfs_cnt = read_prnt_directory(&an_dir[i], i, &sn);
        validate_VL_an_dir_contiguity(&an_dir[i]);
    }
    
    //validate the compatibility of the analysis directories
    validate_VL_an_dir_compatibility(&an_dir[0], dir_count);
    
    //parse input sample names and generate merged sample name
    //an auto-generated sample name is always constructed because
    //the compatibility of sample metadata is validated during
    //this process
    generate_VL_sample_name(&sn);
    
    int mod_detected = 0; //flag that modified  reads were detected
    int unt_detected = 0; //flag that untreated reads were detected
    int den_detected = 0; //flag that denatured reads were detected
        
    //store profiles and normalize reactivity values using whole dataset
    for (i = 0; i < dir_count; i++) {   //for each input analysis directory
        
        mod_detected = unt_detected = den_detected = 0; //zero channel detection variables
        
        for (j = an_dir[i].min_tl; j <= an_dir[i].max_tl; j++) { //for each transcript length
            
            if (an_dir[i].opnd[j]) { //if a profile was opened for the current transcript length
                
                //store the shapemapper2 profile and add the number of
                //target RNA nucleotide reactivity values to the total
                store_SM2_profile(&an_dir[i].data[j], an_dir[i].loc[j]);
                an_dir[i].trg_rct_cnt += an_dir[i].data[j].trg_nt_cnt;
                
                if (an_dir[i].data[j].chnls.mod) {
                    mod_detected = 1;
                }
                
                if (an_dir[i].data[j].chnls.unt) {
                    unt_detected = 1;
                }
                
                if (an_dir[i].data[j].chnls.den) {
                    den_detected = 1;
                }
            }
        }
        
        //set minimum and maximum opened profiles
        set_opnd_profile_bounds(&an_dir[i]);
        
        //set the channel configuration of the entire dataset, confirm that
        //the channel configuration is valid, and confirm that all input
        //datasets have the same channel configuration
        an_dir[i].chnls.mod = mod_detected;
        an_dir[i].chnls.unt = unt_detected;
        an_dir[i].chnls.den = den_detected;
        validate_channel_configuration(&an_dir[i].chnls);
        validate_channel_compatibility(&an_dir[i].chnls, &an_dir[0].chnls);
        
        //set the start index of the entire dataset and confirm that all
        //individual transcripts of the dataset have the same start index
        //and that all analysis directories have the same start index
        an_dir[i].trgt_start = validate_trgt_start(&an_dir[i]);
        validate_ext_start_ix_compatibility(an_dir[i].trgt_start, an_dir[0].trgt_start);
        
        //verify that each transcript sequence is a substring of the
        //next transcript sequence
        validate_transcript_substrings(&an_dir[i]);
        
        //initialize empty profiles for missing transcript lengths
        for (j = an_dir[i].min_tl; j <= an_dir[i].max_tl; j++) {
            if (!an_dir[i].opnd[j]) {
                allocate_SM2_profile_memory(&an_dir[i].data[j], j+an_dir[i].trgt_start);
                initialize_empty_profile(&an_dir[i].data[j], j, an_dir[i].trgt_start);
            }
        }
        
        //perform whole dataset reactivity normalization
        normalize_VL_reactivities(&an_dir[i], min_depth, max_bkg, verify_norm);
    }
    
    //print sample names to screen
    printf("\nuser-specified sample name: %s\n", (sn.usr[0]) ? sn.usr : "not provided");
    printf("auto-generated sample name: %s\n\n", sn.mrg);

    make_VL_output_directories(&an_dir[0], &outfiles, &sn);  //make output directories/files
    
    SM2_analysis_directory mrg = {{0}};          //storage for merged replicate data
    SM2_analysis_directory * data2output = NULL; //pointer to data set to use for output
    
    if (dir_count == 1) {     //if there is only one input directory
        normalize_VL_reactivities(&an_dir[0], min_depth, max_bkg, 0); //normalize input data
        data2output = &an_dir[0];                                     //output that input directory
        
    } else { //if there is more than one input directory
        merge_VL_profiles(&an_dir[0], dir_count, &mrg, min_depth, max_bkg); //merge SM2 profiles
        normalize_VL_reactivities(&mrg, min_depth, max_bkg, 0);             //normalize merged data
        data2output = &mrg;
    }
    
    print_processing_record(&sn, &outfiles, &an_dir[0], &mrg); //print processing record
    print_merged_VL_profiles(data2output, &outfiles);          //print profile output
    print_legacy_compiled_table(data2output, &outfiles, &sn);  //print legacy compiled table
    
    //close output files
    for (i = an_dir[0].min_tl; i <= an_dir[0].max_tl; i++) {
        if (fclose(outfiles.ofp[i]) == EOF) {
            printf("main: error - failed to close output file. Aborting program...\n");
            abort();
        }
    }
    
    //close input files
    for (i = 0; i < dir_count; i++) {
        for (j = an_dir[0].min_tl; j <= an_dir[0].max_tl; j++) {
            if (an_dir[i].opnd[j]) {
                if (fclose(an_dir[i].prf[j]) == EOF) {
                    printf("read_SM2out_directory: error - failed to close input file. Aborting program...\n");
                    abort();
                }
            }
        }
    }
}

/* print_processing_record: print record of input files, merged sample
 name,and transcrip length inventory to screen and file*/
void print_processing_record(sample_names * sn, output_files * outfiles, SM2_analysis_directory * an_dir, SM2_analysis_directory * mrg)
{
    int i = 0; //general purpose index
    
    FILE * p_prcs_rcrd = NULL;         //processing record file pointer
    char prcs_rcrd_nm[MAX_LINE] = {0}; //processing record file name
    char tmp_sn[MAX_LINE] = {0};       //temp sample name
    
    //generate processing record file name
    if (snprintf(prcs_rcrd_nm, MAX_LINE, "./%s/000_processing_record.txt", outfiles->out_dir) >= MAX_LINE) {
        printf("print_processing_record: error - processing record file name exceeded buffer. aborting...\n");
        abort();
    }
    
    //open processing record file
    if ((p_prcs_rcrd = fopen(prcs_rcrd_nm, "w")) == NULL) {
        printf("print_processing_record: error - could not open processing record output file. Aborting program...\n");
        abort();
    }
    
    
    //print input directory information
    printf("input directories:\n");
    fprintf(p_prcs_rcrd, "input directories\n");
    
    for (i = 0; i < sn->cnt; i++) {
        strcpy(tmp_sn, sn->ipt[i]); //copy iput file name to tmp_sn array
        remove_out_suffix(tmp_sn);  //remove transcript length out suffix
                
        //print input sample name to file and screen
        printf("input %d: %s\n", i+1, tmp_sn);
        fprintf(p_prcs_rcrd, "input %d: %s\n", i+1, tmp_sn);
    }
    
    //if more than one input was provided, print merged sample name
    if (sn->cnt > 1) {
        printf("merged:  %s\n", sn->sn2use);
        fprintf(p_prcs_rcrd, "merged:  %s\n", sn->sn2use);
    }
    
    printf("\n");
    fprintf(p_prcs_rcrd, "\n");
    
    
    //print transcript length inventory
    printf("transcript length inventory:\n");
    fprintf(p_prcs_rcrd, "transcript length inventory:\n");
    for (i = 0; i < sn->cnt; i++) {
        printf("input %d: %3d/%d input directories contained a profile\n", i, an_dir[i].prfs_cnt, an_dir[i].outs_cnt);
        fprintf(p_prcs_rcrd, "input %d: %3d/%d input directories contained a profile\n", i, an_dir[i].prfs_cnt, an_dir[i].outs_cnt);
    }
    
    printf("minimum transcript length: %d\n", an_dir[0].min_tl);
    printf("maximum transcript length: %d\n\n", an_dir[0].max_tl);
    
    fprintf(p_prcs_rcrd, "minimum transcript length: %d\n", an_dir[0].min_tl);
    fprintf(p_prcs_rcrd, "maximum transcript length: %d\n\n", an_dir[0].max_tl);
    
    //print normalization factors
    /*printf("normalization factors:\n");
    fprintf(p_prcs_rcrd, "normalization factors\n");
    
    for (i = 0; i < sn->cnt; i++) {
        printf("input %d: %.10f\n", i+1, an_dir[i].cnf);
        fprintf(p_prcs_rcrd, "input %d: %.10f\n", i+1, an_dir[i].cnf);
    }
    
    if (sn->cnt > 1) {
        printf("merged:  %.10f\n", mrg->cnf);
        fprintf(p_prcs_rcrd, "merged:  %.10f\n", mrg->cnf);
    }*/
    
    //printf("\n");
    //fprintf(p_prcs_rcrd, "\n");
    
    //close processing record file
    if (fclose(p_prcs_rcrd) == EOF) {
        printf("print_processing_record: error - failed to close processing record file. Aborting program...\n");
        abort();
    }
}




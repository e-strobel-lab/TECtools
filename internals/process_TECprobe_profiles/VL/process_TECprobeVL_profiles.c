//
//  process_TECprobeVL_profiles.c
//  
//
//  Created by Eric Strobel on 1/19/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include <dirent.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"
#include "../../utils/io_management.h"

#include "../../cotrans_preprocessor/run_script_gen/UNV/config_struct.h"

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
    extern const char empty_SM2out[6];
    
    int mode = -1; //data processing mode
    
    char * prnt_dir_nm[MAX_RUNS] = {NULL}; //list of parent directories
    
    int dir_count = 0; //flag that input directory was supplied
    
    int min_depth = 5000;  //minimum effective read depth
                           //NOTE: SM2 docs say default 5000, but 1000 is used in code
    double max_bkg = 0.05; //maximum untreated mutation rate
    
    int min_depth_provided = 0; //flag that minimum depth option was provided
    int max_bkg_provided = 0;   //flag that maximum background option was provided
    int norm_all = 0;           //flag to normalize non-HQ nucleotides
    int verify_norm = 0;        //flag to verify normalization against SM2 calculations
    
    sample_names sn = {{{0}}};     //structure for sample name storage/merged name construction
    output_files outfiles = {{0}}; //structure for storing output file pointers and names
    
    char dflt_out_dir_nm[20] = {"dataset_norm_out"};
    
    int * ix = NULL; //pointer for target indices
    
    int ret = 0; //variable for storing snprintf return value
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"mode",           required_argument, 0, 'm'}, //set mode
            {"input",          required_argument, 0, 'i'}, //shapemapper analysis directory
            {"out_dir_name",   required_argument, 0, 'o'}, //output directory name
            {"sample-name",    required_argument, 0, 'n'}, //user-supplied sample name
            {"debug",          no_argument,       0, 'd'}, //debug flag
            {"min-depth",      required_argument, 0, 'e'}, //min effective depth
            {"max-background", required_argument, 0, 'b'}, //max background
            {"norm-all",       no_argument,       0, 'a'}, //normalize non-HQ nucleotides
            {"verify_norm",    no_argument,       0, 'y'}, //verify normalization
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "m:i:o:n:de:b:ay", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
            
            case 'm': //set mode
                if (!strcmp(argv[optind-1], "MULTILENGTH")) {
                    mode = MLT;
                } else if (!strcmp(argv[optind-1], "MULTIPLEX")) {
                    mode = MUX;
                } else {
                    printf("process_TECprobe_profiles: error - unexpected run mode. valid run modes currently include: MULTILENGTH. aborting...\n");
                    abort();
                }
                
                break;
            
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
                ret = snprintf(outfiles.out_dir, MAX_NAME, "%s_%s", argv[optind-1], dflt_out_dir_nm);
                if (ret >= MAX_NAME || ret < 0) {
                    printf("process_TECprobeVL_profiles: error - error when storing output directory name. aborting...\n");
                    abort();
                }
                
                break;
                
            case 'n': //user-supplied output file name
                
                //store user-supplied output file name
                ret = snprintf(sn.usr, MAX_NAME, "%s", argv[optind-1]);
                if (ret >= MAX_NAME || ret < 0) {
                    printf("process_TECprobeVL_profiles: error - error when storing sample name. aborting...\n");
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
                
            case 'a': //normalize non-HQ nucleotides too
                norm_all = 1;
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
    
    if (mode == -1) {
        printf("process_TECprobe_profiles: error - run mode was not set. valid run modes currently include: MULTILENGTH. aborting...\n");
        abort();
    }
    
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
    if (sn.usr[0]) {               //if a user-supplied output sample name was provided
        sn.sn2use = &sn.usr[0];    //use the user-supplied output sample name
        
    } else if (dir_count == 1) {   //otherwise, if only 1 input directory was provided
        sn.sn2use = &sn.ipt[0][0]; //use sample name derived from the input
        
    } else if (dir_count > 1) {    //otherwise, if more than 1 input directory was provided
        sn.sn2use = &sn.mrg[0];    //use the auto-generated sample name
    }
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
        
    SM2_analysis_directory * an_dir = NULL; //pointer for SM2 analysis directory memory allocation
    
    //allocate memory for storing SM2 analysis directory information
    if ((an_dir = calloc(dir_count, sizeof(*an_dir))) == NULL) {
        printf("process_TECprobeVL_profiles: error - memory allocation for directory/file pointers failed. aborting...\n");
        abort();
    }
    
    //initialize SM2 analysis directory structs
    for (i = 0; i < dir_count; i++) {                  //for each input directory
        strcpy(an_dir[i].prnt_dir_nm, prnt_dir_nm[i]); //store the parent directory name
        an_dir[i].min_id = INT_MAX;                    //initialize min_tl
        an_dir[i].max_id = INT_MIN;                    //initialize max_tl
    }
        
    //read each SM2 analysis directory and confirm that the number of
    //files and transcript length range is the same for all samples
    for (i = 0; i < dir_count; i++) {
        
        //open all relevant directories and files in the parent directory
        //then validate that transcript length profiles are contiguous
        an_dir[i].prfs_cnt = read_prnt_directory(&an_dir[i], i, &sn);
        
        if (mode == MLT) {
            validate_VL_an_dir_contiguity(&an_dir[i]);
        }
    }
    
    //validate the compatibility of the analysis directories
    validate_an_dir_compatibility(&an_dir[0], dir_count);
    
    //parse input sample names and generate merged sample name
    //an auto-generated sample name is always constructed because
    //the compatibility of sample metadata is validated during
    //this process
    generate_VL_sample_name(&sn);
    
    int mod_detected[MAX_RUNS] = {0}; //flag that modified  reads were detected
    int unt_detected[MAX_RUNS] = {0}; //flag that untreated reads were detected
    int den_detected[MAX_RUNS] = {0}; //flag that denatured reads were detected
    
    int empty_prf_tot_lines = 0; //number of lines to include in an empty profile
    int empty_prf_trg_lines = 0; //number of target lines in an empty profile
    int set_empty_cnts = 0;      //flag that empty target nucleotide counts were set
        
    //store profiles and normalize reactivity values using whole dataset
    for (i = 0; i < dir_count; i++) {   //for each input analysis directory
        
        ix = &an_dir[i].indx[0]; //set pointer for target indices
        for (j = 0; ix[j] <= an_dir[i].max_id; j++) { //for each target
            if (strcmp(an_dir[i].loc[ix[j]], empty_SM2out)) {
                
                //store the shapemapper2 profile and add the number of
                //target RNA nucleotide reactivity values to the total
                store_SM2_profile(&an_dir[i].data[ix[j]], an_dir[i].loc[ix[j]]);
                an_dir[i].trg_rct_cnt += an_dir[i].data[ix[j]].trg_nt_cnt;
                
                if (an_dir[i].data[ix[j]].chnls.mod) {
                    mod_detected[i] = 1;
                }
                
                if (an_dir[i].data[ix[j]].chnls.unt) {
                    unt_detected[i] = 1;
                }
                
                if (an_dir[i].data[ix[j]].chnls.den) {
                    den_detected[i] = 1;
                }
            }
        }
    }
    
    for (i = 0; i < dir_count; i++) {
        
        for (j = 0; ix[j] <= an_dir[i].max_id; j++) { //for each target
            if (!strcmp(an_dir[i].loc[ix[j]], empty_SM2out)) {
                
                //initialize empty profiles for missing targets
                //target start from the first target in the current directory is used during initialization
                
                empty_prf_tot_lines = 0;
                empty_prf_trg_lines = 0;
                
                if (mode == MLT) {
                    //in multilength mode (TECprobe-VL and TECprobe-LM experiments),
                    //the target id is the target sequence length and can be used to
                    //set the total and target line counts for empty profiles
                    
                    empty_prf_tot_lines = ix[j] + an_dir[i].data[ix[0]].trgt_start;
                    empty_prf_trg_lines = ix[j];
                    
                } else  if (mode == MUX) {
                    //in multiplex mode (TECprobe-MUX experiments), the target id is
                    //the barcode id, and is therefore unrelated to target sequence
                    //length. therefore, the code below attempts to determine total
                    //and target line counts from the corresponding target in a
                    //replicate data set.
                    
                    for (k = 0, set_empty_cnts = 0; k < dir_count && !set_empty_cnts; k++) {
                        if (strcmp(an_dir[k].loc[ix[j]], empty_SM2out)) {
                            empty_prf_tot_lines = an_dir[k].data[ix[j]].tot_nt_cnt;
                            empty_prf_trg_lines = an_dir[k].data[ix[j]].trg_nt_cnt;
                            set_empty_cnts = 1;
                        }
                    }
                    
                    //if no data sets contain a profile for the target, the total line
                    //count is set to the number of nucleotides that precede the target
                    //sequence in other profiles, and the number of target lines is set
                    //to zero.
                    
                    if (!set_empty_cnts) {
                        empty_prf_tot_lines = an_dir[i].data[ix[0]].trgt_start;
                        empty_prf_trg_lines = 0;
                    }
                    
                }
                allocate_SM2_profile_memory(&an_dir[i].data[ix[j]], empty_prf_tot_lines);
                initialize_empty_profile(&an_dir[i].data[ix[j]], empty_prf_trg_lines, an_dir[i].data[ix[0]].trgt_start);
            }
        }
        
        //set minimum and maximum target lengths
        set_len_range(&an_dir[i]);
        
        //set the channel configuration of the entire dataset, confirm that
        //the channel configuration is valid, and confirm that all input
        //datasets have the same channel configuration
        an_dir[i].chnls.mod = mod_detected[i];
        an_dir[i].chnls.unt = unt_detected[i];
        an_dir[i].chnls.den = den_detected[i];
        validate_channel_configuration(&an_dir[i].chnls);
        validate_channel_compatibility(&an_dir[i].chnls, &an_dir[0].chnls);
        
        //set the start index of the entire dataset and confirm that all
        //individual targets of the dataset have the same start index
        //and that all analysis directories have the same start index
        an_dir[i].trgt_start = validate_trgt_start(&an_dir[i]);
        validate_ext_start_ix_compatibility(an_dir[i].trgt_start, an_dir[0].trgt_start);
        
        //verify that each transcript sequence is a substring of the
        //next transcript sequence
        validate_transcript_substrings(&an_dir[i]);
        
        //perform whole dataset reactivity normalization
        normalize_VL_reactivities(&an_dir[i], min_depth, max_bkg, norm_all, verify_norm);
    }
    
    //print sample names to screen
    printf("\nuser-specified sample name: %s\n", (sn.usr[0]) ? sn.usr : "not provided");
    printf("auto-generated sample name: %s\n\n", sn.mrg);
    
    SM2_analysis_directory mrg = {{0}};          //storage for merged replicate data
    SM2_analysis_directory * data2output = NULL; //pointer to data set to use for output
    
    if (dir_count == 1) {     //if there is only one input directory
        normalize_VL_reactivities(&an_dir[0], min_depth, max_bkg, norm_all, 0); //normalize input data
        data2output = &an_dir[0];                                               //output that input directory
        
    } else { //if there is more than one input directory
        merge_VL_profiles(&an_dir[0], dir_count, &mrg, min_depth, max_bkg);     //merge SM2 profiles
        normalize_VL_reactivities(&mrg, min_depth, max_bkg, norm_all, 0);       //normalize merged data
        data2output = &mrg;
    }
    
    if ((outfiles.ofp = calloc(data2output->sd_cnt, sizeof(*outfiles.ofp))) == NULL) {
        printf("process_TECprobe_VL_profiles: error - failed to allocate memory for output files. aborting...\n");
        abort();
    }
    
    make_VL_output_directories(&an_dir[0], &outfiles, &sn);  //make output directories/files
    
    print_processing_record(&sn, &outfiles, &an_dir[0], &mrg); //print processing record
    print_merged_VL_profiles(data2output, &outfiles);          //print profile output
    print_legacy_compiled_table(data2output, &outfiles, &sn);  //print legacy compiled table
    
    //close output files
    ix = &data2output->indx[0]; //set pointer to target indices
    for (i = 0; ix[i] <= data2output->max_id; i++) {
        if (fclose(outfiles.ofp[ix[i]]) == EOF) {
            printf("main: error - failed to close output file. Aborting program...\n");
            abort();
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
    
    int ret = 0; //variable for storing snprintf return value
    
    //generate processing record file name
    ret = snprintf(prcs_rcrd_nm, MAX_LINE, "./%s/000_processing_record.txt", outfiles->out_dir);
    if (ret >= MAX_LINE || ret < 0) {
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
        //print input sample name to file and screen
        printf("input %d: %s\n", i+1, sn->ipt[i]);
        fprintf(p_prcs_rcrd, "input %d: %s\n", i+1, sn->ipt[i]);
    }
    
    //if more than one input was provided, print merged sample name
    if (sn->cnt > 1) {
        printf("merged:  %s\n", sn->sn2use);
        fprintf(p_prcs_rcrd, "merged:  %s\n", sn->sn2use);
    }
    
    printf("\n");
    fprintf(p_prcs_rcrd, "\n");
    
    
    //print target inventory
    printf("target inventory:\n");
    fprintf(p_prcs_rcrd, "target inventory:\n");
    for (i = 0; i < sn->cnt; i++) {
        printf("input %d: %3d/%d input directories contained a profile\n", i+1, an_dir[i].prfs_cnt, an_dir[i].outs_cnt);
        fprintf(p_prcs_rcrd, "input %d: %3d/%d input directories contained a profile\n", i+1, an_dir[i].prfs_cnt, an_dir[i].outs_cnt);
    }
    
    printf("minimum transcript length: %d\n", an_dir[0].len[MIN]);
    printf("maximum transcript length: %d\n\n", an_dir[0].len[MAX]);
    
    fprintf(p_prcs_rcrd, "minimum transcript length: %d\n", an_dir[0].len[MIN]);
    fprintf(p_prcs_rcrd, "maximum transcript length: %d\n\n", an_dir[0].len[MAX]);
    
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




//
//  assemble_TECprobeLM_data.c
//  
//
//  Created by Eric Strobel on 7/5/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <string.h>

#include "../utils/io_management.h"
#include "../global/global_defs.h"
#include "../mkmtrx/cotrans_mtrx.h"

#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

#include "./store_ipt_name.h"
#include "./validate_input.h"
#include "./count_delims_2_col.h"
#include "./get_value.h"

#include "./print_input_filenames.h"
#include "./print_reactivity_output.h"
#include "./print_length_dist_output.h"

#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"


int main(int argc, char *argv[])
{
    input_data ipt = {{{0}}};
    
    mode_parameters mode_params = {0};             //mode parameters
    mode_params.mod = MODE_INIT;                   //intialize mode
    
    int nrchd_len[TOT_SAMPLES] = {0};              //array to store enriched transcript length values (one per sample)
    int delims2col[TOT_SAMPLES][MAX_IPT] = {0};    //number of delimiters to reach the target data column
    
    double vals[MAX_TRANSCRIPT][TOT_SAMPLES][MAX_IPT] = {{{0}}}; //array to store values from the target data columns
    
    char out_nm[MAX_NAME] = {0};    //output file name
    char out_dir[MAX_NAME*2] = {0}; //output directory name
        
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    
    while (1) {
        static struct option long_options[] =
        {
            {"mode",       required_argument,  0,  'm'}, //run mode (REACTIVITY or LEN_DIST)
            {"enriched1",  required_argument,  0,  'a'}, //enriched length for sample 1
            {"enriched2",  required_argument,  0,  'b'}, //enriched length for sample 2
            {"enriched3",  required_argument,  0,  'c'}, //enriched length for sample 3
            {"data1",      required_argument,  0,  'x'}, //input file for sample 1
            {"data2",      required_argument,  0,  'y'}, //input file for sample 2
            {"data3",      required_argument,  0,  'z'}, //input file for sample 3
            {"out_name",   required_argument,  0,  'o'}, //output file name
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "m:a:b:c:x:y:z:o:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
              
            case 'm': //set mode
                if (mode_params.mod != -1) { //mode was already set, throw error
                    printf("assemble_TECprobeLM_data: error - more than one mode argument was provided. aborting...\n");
                    abort();
                } else {
                    if (!strcmp(argv[optind-1], "REACTIVITY")) {                  //REACTIVITY mode
                        mode_params.mod = REACTIVITY;                             //set mode to REACTIVITY
                        mode_params.dlm = ',';                                    //set delimiter to comma
                        mode_params.offset = 1;                                   //set offset to 1
                        strcpy(mode_params.hdr, "Norm_calc_profile_");            //set target header to "Norm_calc_profile_"
                    } else if (!strcmp(argv[optind-1], "LEN_DIST")) {  //LEN_DIST mode
                        mode_params.mod = LEN_DIST;                    //set mode to LEN_DIST
                        mode_params.dlm = '\t';                                   //set delimiter to tab
                        //NOTE: offset is set later based on the minimum transcript length
                        strcpy(mode_params.hdr, "frac_algnd");                    //set target header to frac_algnd
                    } else {
                        printf("assemble_TECprobeLM_data: error - unrecognized mode. aborting...\n");
                        abort();
                    }
                }
                break;
                
            case 'a': //set sample 1 enriched length
                if (nrchd_len[S1]) { //enriched length 1 was already set, throw error
                    printf("assemble_TECprobeLM_data: error - more than one value provided for sample 1 enriched length. aborting...\n");
                    abort();
                } else {
                    nrchd_len[S1] = atoi(argv[optind-1]);
                }
                
                break;
                
            case 'b': //set sample 2 enriched length
                if (nrchd_len[S2]) { //enriched length 2 was already set, throw error
                    printf("assemble_TECprobeLM_data: error - more than one value provided for sample 2 enriched length. aborting...\n");
                    abort();
                } else {
                    nrchd_len[S2] = atoi(argv[optind-1]);
                }
                
                break;
                
            case 'c': //set sample 3 enriched length
                if (nrchd_len[S3]) { //enriched length 3 was already set, throw error
                    printf("assemble_TECprobeLM_data: error - more than one value provided for sample 3 enriched length. aborting...\n");
                    abort();
                } else {
                    nrchd_len[S3] = atoi(argv[optind-1]);
                }
                break;
                
            case 'x': //set sample 1 input
                if (ipt.cnt[S1] < MAX_IPT) { //check that max input hasn't been reached
                    get_file(&(ipt.fp[S1][ipt.cnt[S1]]), argv[optind-1]);
                    store_ipt_name(&ipt.fn[S1][ipt.cnt[S1]], argv[optind-1]);
                    ipt.cnt[S1]++;
                } else {
                    printf("assemble_TECprobeLM_data: error too many sample 1 input files provided. aborting...");
                    abort();
                }
                break;
            
            case 'y': //set sample 2 input
                if (ipt.cnt[S2] < MAX_IPT) { //check that max input hasn't been reached
                    get_file(&(ipt.fp[S2][ipt.cnt[S2]]), argv[optind-1]);
                    store_ipt_name(&ipt.fn[S2][ipt.cnt[S2]], argv[optind-1]);
                    ipt.cnt[S2]++;
                } else {
                    printf("assemble_TECprobeLM_data: error too many sample 2 input files provided. aborting...");
                    abort();
                }
                break;
            
            case 'z': //set sample 3 input
                if (ipt.cnt[S3] < MAX_IPT) { //check that max input hasn't been reached
                    get_file(&(ipt.fp[S3][ipt.cnt[S3]]), argv[optind-1]);
                    store_ipt_name(&ipt.fn[S3][ipt.cnt[S3]], argv[optind-1]);
                    ipt.cnt[S3]++;
                } else {
                    printf("assemble_TECprobeLM_data: error too many sample 3 input files provided. aborting...");
                    abort();
                }
                break;
                
            case 'o': //set output file name
                if (!out_nm[0]) { //check that output name hasn't been set
                    strcpy(out_nm, argv[optind-1]);
                } else {
                    printf("assemble_TECprobeLM_data: error more than one output name provided. aborting...");
                    abort();
                }
                
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
    
    //check that mode was set
    if (mode_params.mod == -1) {
        printf("assemble_TECprobeLM_data: error - run mode was not set. set run mode to REACTIVITY or LEN_DIST depending on input file type. aborting...\n");
        abort();
    }
    
    //validate input files
    validate_input(&ipt, mode_params.mod);
    
    //check that output name was provided and make output directory
    if (!out_nm[0]) {
        printf("assemble_TECprobeLM_data: error - no output name was provided. aborting...\n");
        abort();
    } else {
        if (mode_params.mod == REACTIVITY) {
            sprintf(out_dir, "%s_reactivity", out_nm);
        } else if (mode_params.mod == LEN_DIST) {
            sprintf(out_dir, "%s_lenDist", out_nm);
        } else {
            abort();
        }
        mk_out_dir(out_dir);
    }
    
    //print record of input file names
    print_input_filenames(out_dir, out_nm, ipt.fn, ipt.cnt);

    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
        
    //for every input file of every sample, count the number of delimiters
    //it takes to reach the column that contains the desired values
    for (i = 0; i < TOT_SAMPLES; i++) {
        for (j = 0; j < ipt.cnt[i]; j++) {
            if (!count_delims_2_col(ipt.fp[i][j], &mode_params, nrchd_len[i], &(delims2col[i][j]))) {
                printf("assemble_TECprobeLM_data: error - delim count failed. aborting...\n");
                abort();
            }
        }
    }
    
    int found_end = 0; //flag that end of file was found
    int line_term = 0; //terminating character of the current line
    int ipt0_term = 0; //terminating character of the current line in the first input file that was assessed
    int max_index = 0; //maximum index after parsing alignment rate file (== max transcript length - 1)
    
    for (i = 0; i < MAX_TRANSCRIPT && (i < nrchd_len[S3] || mode_params.mod == LEN_DIST) && !found_end; i++) {
        //the loop is run until either:
        //1. MAX_TRANSCRIPT is reached
        //2. the maximum enriched length of sample 3 is reach (REACTIVITY mode only)
        //3. the end of the file is found
        
        for (j = 0; j < TOT_SAMPLES; j++) { //for each sample
            
            if (i < nrchd_len[j] || mode_params.mod == LEN_DIST) {
                //perform the operations below if either:
                //1. i is less than the enriched length for the current sample (REACTIVITY mode only)
                //2. LEN_DIST mode is being run
                
                //note: the reactivity of the last nt is always zero,
                //but this zero is not in compiled matrix file
                
                for (k = 0; k < ipt.cnt[j]; k++) { //for each input file of the current sample
                    
                    //get the value of the target column from the current line
                    if ((line_term = get_value(ipt.fp[j][k], &mode_params, delims2col[j][k], &(vals[i][j][k])))) {
                        
                        if (k == 0) {               //if processing the first input file of the current sample
                            ipt0_term = line_term;  //set the terminating character of the current line
                            if (line_term == EOF) { //if the terminating character is EOF
                                found_end = 1;      //set flag that the end of the file was found
                            }
                        } else {                    //if processing a non-first input file
                            if (line_term != ipt0_term) { //check that the terminating character matches that of the first input file
                                printf("assemble_TECprobeLM_data: error - lines at same file position terminate with different characters\n");
                                abort();
                            }
                        }
                        
                        if (!found_end) {                             //if the end of the file was not found
                            if (j == 0 && k == 0) {                   //if processing the first input file of the first sample
                                printf("%d", i+mode_params.offset); //print the id of the current line (nucleotide or transcript)
                            }
                            printf("\t%10.6f", vals[i][j][k]);        //print the value of the target column
                        }
                        
                    } else { //throw error if no value was found
                        printf("\n\nerror: no value found. aborting...\n");
                        abort();
                    }
                }
                
            } else { //if there is no target value to get
                
                if (j == 0) {                             //if processing the first sample
                    printf("%d", i+mode_params.offset);   //print the id of the current line (nucleotide or transcript
                }
                
                for (k = 0; k < ipt.cnt[j]; k++) {        //for every input file of the current sample
                    printf("\t    --    ");               //print a spacer
                }
            }
        }
        printf("\n");
    }
    
    if (mode_params.mod == REACTIVITY) { //if in reactivity mode, print reactivity output
        print_reactivity_output(out_dir, out_nm, &mode_params, ipt.cnt, nrchd_len, vals);
        
    } else if (mode_params.mod == LEN_DIST) { //if in LEN_DIST mode, print length distribution output
        max_index = i - 1; //set maximum data-containing line index
        print_length_dist_output(out_dir, out_nm, &mode_params, ipt.cnt, max_index, vals);
    }
    
}

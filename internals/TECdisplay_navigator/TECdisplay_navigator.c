//
//  TECdisplay_navigator.c
//  
//
//  Created by Eric Strobel on 7/2/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <sys/stat.h>
#include <math.h>

#include "../utils/io_management.h"
#include "../utils/debug.h"

#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/ispair.h"
#include "../seq_utils/is_dgnrt_mtch.h"
#include "../seq_utils/test_possible_pairs.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"
#include "./merge_values_files.h"
#include "./parse_reference.h"
#include "./parse_constraints.h"
#include "./filter_values.h"

int main(int argc, char *argv[])
{
    extern int debug; //flag to run debug mode
    
    values_input vals[MAX_VALS] = {0};         //storage for sample name, file name, and file pointer of values files
    FILE * fp_cons = NULL;                     //pointer to constraints file
    
    int vals_cnt = 0;                          //number of values files provided
    int cons_provided = 0;                     //number of constraint files provided
    int out_nm_provided = 0;                   //number of output directory names provided
    
    char out_dir_nm[MAX_NAME] = {0};           //array to store output directory name
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"values",        required_argument,  0,  'v'},  //values file input
            {"constraints",   required_argument,  0,  'c'},  //constraints file input
            {"debug",   	  required_argument,  0,  'd'},  //turn on debug mode
            {"out-name",      required_argument,  0,  'o'},  //set output directory name
            
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "v:c:d:o:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
           
            case 'v': //get values file
                if (vals_cnt < MAX_VALS) {                              //check that MAX_VALS will not be exceeded
                    strcpy(vals[vals_cnt].fn, argv[optind-1]);          //store file name
                    get_sample_name(argv[optind-1], vals[vals_cnt].nm); //get sample name from file name
                    get_file(&(vals[vals_cnt].fp), argv[optind-1]);     //set file pointer to values file
                    vals_cnt++;                                         //count values files provided
                } else {
                    printf("TECdisplay_navigator: error - too many values files were provided. the maximum is %d. aborting...\n", MAX_VALS);
                    abort();
                }
                break;
            
            case 'c': //get constraints file
                if (!cons_provided) {                     //check that a constraints file was not previously provided
                    get_file(&(fp_cons), argv[optind-1]); //set file pointer to constraints file
                    cons_provided++;                      //count constraints files provided
                } else {
                    printf("TECdisplay_navigator: error - more than one constraints file was provided. aborting..\n");
                    abort();
                }
                break;
                
            case 'o': //set output_directory name
                if (!out_nm_provided) {                        //check that output directory was not previously provided
                    if (strlen(argv[optind-1]) < (MAX_NAME)) { //check output directory name length
                        strcpy(out_dir_nm, argv[optind-1]);    //store output directory name
                    } else {
                        printf("TECdisplay_navigator: error - output directory name is longer than the maximum length (%d). setting output directory to default name 'out'.\n", MAX_NAME);
                    }
                    out_nm_provided++; //increment output name counter
                    
                } else {
                    printf("TECdisplay_navigator: error - more than one output directory name was provided. aborting\n");
                    abort();
                }
                break;
            
            case 'd': //turn on debug mode
                debug = 1;
                break;
                
            default:
                printf("error: unrecognized option. Aborting program...\n");
                abort();
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

    /* make output directory */
    if (!out_dir_nm[0]) {          //if no output directory name was provided
        strcpy(out_dir_nm, "out"); //if not, set output directory name to "out"
    }
    mk_out_dir(out_dir_nm);        //make output directory
    
    
    /* merge multiple input files into a single file. if only one input file is supplied, it is still
     processed by the merge_values_files function to append the sample name to the column headers. as
     written, this is required in the get_sample_info function later */
    char merged_out_nm[MAX_LINE] = {0}; //merged output filename
    
    if (!out_dir_nm[0]) {                                        //if no output directory name was provided
        sprintf(merged_out_nm, "%s/merged_out.txt", out_dir_nm); //set merged values file name to merged_out.txt
    } else {                                                                    //otherwise
        sprintf(merged_out_nm, "%s/%s_merged_out.txt", out_dir_nm, out_dir_nm); //set to <out_dir_nm>_merged_out.txt
    }
    merge_values_files(&vals[0], vals_cnt, merged_out_nm);       //merge input files
    
    FILE * ipt = NULL;                               //pointer to merged values file
    if ((ipt = fopen(merged_out_nm, "r")) == NULL) { //open merged values file
        printf("TECdisplay_navigator: error - could not open merged input file. aborting...\n");
        abort();
    }
    
    wt_source wt = {0};                        //wt source sequence
    basemap bmap = {0};                        //basemap for storing variable base reference sequence
    char * cnstnt_indels = NULL;               //storage for constant indels string
    constraints cons[MAX_CONSTRAINTS] = {{0}}; //constraint parameters
    int cons_cnt = 0;                          //number of constraints in constraints file
    
    parse_reference(fp_cons, &bmap, &wt, &cnstnt_indels);                  //construct basemap from reference sequence
    cons_cnt = parse_constraints(fp_cons, &cons[0], &bmap, cnstnt_indels); //parse constraints file and set constraints
    filter_values(ipt, &cons[0], cons_cnt, &bmap, out_dir_nm);             //filter values for matches to constraints
    
    /* close merged values file */
    if (fclose(ipt) == EOF) {
        printf("TECdisplay_navigator: error - error occurred when closing merged input file. Aborting program...\n");
        abort();
    }
    
    return 1;
}















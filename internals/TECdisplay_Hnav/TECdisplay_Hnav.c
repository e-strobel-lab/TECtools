//
//  TECdisplay_Hnav.c
//  
//
//  Created by Eric Strobel on 11/3/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/stat.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"
#include "../utils/debug.h"

#include "../seq_utils/basemap.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/test_possible_pairs.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"

#include "../TECdisplay_navigator/TECdisplay_navigator_defs.h"
#include "../TECdisplay_navigator/TECdisplay_navigator_structs.h"
#include "../TECdisplay_navigator/merge_values_files.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

#include "./get_constraint_metadata.h"
#include "./calc_output_files.h"
#include "./Hfilter.h"
#include "process_output_files.h"

int main(int argc, char *argv[])
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    extern char TDHN_TECDnav_path[MAX_LINE];      //TECdisplay_navigator path
    extern char TDHN_merged_out_nm[15];           //merged_output filename
    extern output_file_names out_fns[MAX_LAYERS]; //storage for TECdisplay_navigator output filenames
    
    extern int debug; //flag to run debug mode
    
    values_input vals[MAX_VALS] = {0};         //storage for sample name, file name, and file pointer of values files
    
    int vals_cnt = 0;                          //number of values files provided
    int layr_cnt = 0;                          //number of constraint file layers provided
    int out_nm_provided = 0;                   //number of output directory names provided
    int out_prfx_provided = 0;                 //number of output prefixes provided
    int nonstandard = 1;                       //flag that input is non-standard format
    int exclude = 0;                           //flag to exclude variants that match any constraint
    char col_id[MAX_NAME] = {0};               //column to print to aggregated output file TODO: add option to set this
    
    char out_dir_nm[MAX_NAME] = {0};           //array to store output directory name
    char out_prefix[MAX_NAME] = {0};           //array to store output prefix
        
    constraint_metadata cons_meta[MAX_LAYERS] = {0}; //storage for constraint file metadata
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"values",        required_argument,  0,  'v'},  //values file input
            {"constraints",   required_argument,  0,  'c'},  //constraints file input
            {"exclude",       required_argument,  0,  'x'},  //exclusion constraints file input
            {"debug",         required_argument,  0,  'd'},  //turn on debug mode
            {"out-name",      required_argument,  0,  'o'},  //set output directory name
            {"out-prefix",    required_argument,  0,  'f'},  //set output file prefix
            {"path",          required_argument,  0,  'p'},  //set TECdisplay_navigator path
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "v:c:x:d:o:f:p:", long_options, &option_index);
        
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
                if (layr_cnt < MAX_LAYERS) {                                      //check that MAX_LAYERS was not exceeded
                    get_constraint_metadata(argv[optind-1], &cons_meta[layr_cnt], 'c', layr_cnt); //process the contraints file
                    layr_cnt++;                                                         //count constraints files provided
                } else {
                    printf("TECdisplay_navigator: error - too many constraints files were provied. the maximum number of layers is %d. aborting..\n", MAX_LAYERS);
                    abort();
                }
                break;
                
            case 'x': //get exclusion constraints file
                if (layr_cnt < MAX_LAYERS) {                                      //check that MAX_LAYERS was not exceeded
                    get_constraint_metadata(argv[optind-1], &cons_meta[layr_cnt], 'x', layr_cnt); //process the contraints file
                    layr_cnt++;                                                         //count constraints files provided
                } else {
                    printf("TECdisplay_navigator: error - too many constraints files were provied. the maximum number of layers is %d. aborting..\n", MAX_LAYERS);
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
                
            case 'f': //set output file prefix
                if (!out_prfx_provided) {                      //check that output prefix was not previously provided
                    if (strlen(argv[optind-1]) < (MAX_NAME)) { //check output prefix name length
                        strcpy(out_prefix, argv[optind-1]);    //store output prefix name
                    } else {
                        printf("TECdisplay_navigator: error - output prefix name is longer than the maximum length (%d). aborting...\n", MAX_NAME);
                        abort();
                    }
                    out_prfx_provided++; //increment output name counter
                    
                } else {
                    printf("TECdisplay_navigator: error - more than one output prefix name was provided. aborting\n");
                    abort();
                }
                break;
                
            case 'd': //turn on debug mode
                debug = 1;
                break;
                
            case 'p': //set TECdisplay_navigator path
                strcpy(TDHN_TECDnav_path, argv[optind-1]);
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
    
    
    /**** TODO: ADD CHECK THAT ALL CONSTRAINT FILE HEADERS MATCH ****/
    
    if (!out_prfx_provided) { //output prefix is required
        printf("TECdisplay_Hnav: error - no output file prefix was supplied. please provide an output file prefix that describes the input files using the -f option. aborting...\n");
        abort();
    }
    
    if (!TDHN_TECDnav_path[0]) {                            //if a path to TECdisplay_navigator wasn't provided
        strcpy(TDHN_TECDnav_path, "TECdisplay_navigator");  //set path to "TECdisplay_navigator"
    }
    
    if (!col_id[0]) {                                       //if a column id was not provided
        strcpy(col_id, TECdsply_clmn_hdrs[TDSPLY_FRC_HDR]); //set the column to aggregate as "fracBound"
    }
    
    /* make output directory */
    if (!out_dir_nm[0]) {          //if no output directory name was provided
        strcpy(out_dir_nm, "out"); //set output directory name to "out"
    }
    mk_out_dir(out_dir_nm);        //make output directory
    chdir(out_dir_nm);             //change to output directory
    
    print_constraint_metadata(layr_cnt, &cons_meta[0]);
    valid8_constraint_compatiblity(layr_cnt, &cons_meta[0]);
    
    char * layr_list[MAX_LAYERS] = {NULL}; //array to store layer names for temp output directory name construction
    
    calc_output_files(layr_cnt, &cons_meta[0]);                                        //calculate the number of output files
    merge_values_files(&vals[0], vals_cnt, TDHN_merged_out_nm, nonstandard);           //merge input values files into 1 file
    Hfilter(&cons_meta[0], out_prefix, 0, layr_cnt, NULL, layr_list);                  //hierarchically filter the data
    clean_up_output(layr_cnt, &cons_meta[0]);                                          //clean up output directories
    aggregate_output(vals_cnt, &vals[0], layr_cnt, &cons_meta[0], out_prefix, col_id); //make aggregate output files
    
    //print list of TECdisplay_navigator output files for each layer
    int i = 0; //general purpose index
    int j = 0; //general purpose indexs
    
    for (i = 0; i < layr_cnt; i++) {
        printf("layer %d output files:\n", i);
        for (j = 0; j < out_fns[i].f_cnt; j++) {
            printf("%s\n", out_fns[i].fn[j]);
        }
        printf("\n");
    }
}


//
//  mkmtrx.c
//  
//
//  Created by Eric Strobel on 9/24/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include <dirent.h>

#include "../global/global_defs.h"
#include "../utils/io_management.h"
#include "../utils/gen_utils.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"

#include "./get_smpl_nm.h"
#include "./get_fastp_out.h"
#include "./get_profiles.h"
#include "./parse_log.h"
#include "./mk_output_files.h"

#include "mk_rdat.h"

int main(int argc, char *argv[])
{
    FILE * fp_config = NULL;
    
    char prnt_dir[MAX_LINE] = {0};       //parent directory
    
    cotrans_matrix mtrx = {{0}};           //cotrans_matrix struct
    alignment_stats algn[MAX_ROW] = {{{0}}}; //alignment_stats struct
    
    int incld_up2 = 0;          //transcript length to include up to when generating matrix
    int excld_trmnl = 0;        //exclude transcripts at and beyond this position when generating matrix
    char mode_str[16] = {0};    //string to store mode option
    int mode = -1;              //mode value
    
    int mode_set = 0;           //flag that mode was set
    int dir_provided = 0;       //flag that input directory was supplied
    int incUp2_provided = 0;    //flag that include-up-to option was supplied
    int excTrm_provided = 0;    //flag that exclude-term option was supplied
    
    int test_SM2_data = 0;          //flag that turns of fastp output and log parsing when using test SM2 data
    
    int rdat_config_provided = 0; //flag that rdat config was provided
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    
    while (1) {
        static struct option long_options[] =
        {
            {"mode",          required_argument,  0,  'm'},  //shapemapper analysis directory input
            {"input",         required_argument,  0,  'i'},  //shapemapper analysis directory input
            {"include-up-to", required_argument,  0,  'w'},  //whitelist up to arg for transcript length inclusion
            {"exclude-term",  no_argument,        0,  'x'},  //exclude transcripts after arg
            {"test_SM2_data", no_argument,        0,  't'},  //flag that test shapemapper2 data is being assessed
            {"make_rdat",     required_argument,  0,  'r'},  //make rdat file
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "m:i:w:xtr:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
            case 'm': //mode
                if (!mode_set) {
                    strcpy(mode_str, argv[optind-1]); //store input directory name in prnt_dir array
                    if (!strcmp(mode_str, "MULTI")) {
                        mode = MULTI;
                    } else if (!strcmp(mode_str, "SINGLE")) {
                        mode = SINGLE;
                    } else {
                        printf("mkmtrx: error - unrecognized mode\n");
                        abort();
                    }
                    mode_set++;
                } else {
                    printf("mkmtrx: error - more than one mode option provided\n");
                    abort();
                }
                break;
                
            case 'i': //input directory
                if (!dir_provided) {
                    strcpy(prnt_dir, argv[optind-1]); //store input directory name in prnt_dir array
                    dir_provided++;
                } else {
                    printf("mkmtrx: error - more than one input directory was provided. aborting...\n");
                    abort();
                }
                
                break;
                
            case 'w': //whitelist transcripts up to arg
                if (!incld_up2) {
                    incld_up2 = atoi(argv[optind-1]); //set include-up-to limit
                    incUp2_provided++;
                } else {
                    printf("mkmtrx: error - more than one include-up-to was provided. aborting...\n");
                    abort();
                }
                
                break;
                
            case 'x': //exclude transcripts after arg
                if (!excTrm_provided) {
                    excld_trmnl = 1; //set exclude-terminal to true
                    excTrm_provided++;
                } else {
                    printf("mkmtrx: error - more than one exclude-term was provided. aborting...\n");
                    abort();
                }
                
                break;
            
            case 't': //using test SM2 data, disables alignment analyses
                test_SM2_data = 1;
                break;
                
            case 'r': //make rdat file
                if (!rdat_config_provided) {
                    get_file(&(fp_config), argv[optind-1]); //set file pointer to input rdat config file
                    rdat_config_provided++;                 //count rdat config files provided
                } else {
                    printf("mkmtrx: error - more than one rdat config file was provided. aborting...\n");
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

    //check that mode was specified
    if (!mode_set) {
        printf("mkmtrx: error - mode was not set. aborting...\n");
        abort();
    }
    
    //check that input directory was provided
    if (!dir_provided) {
        printf("mkmtrx: error - no input directory was provided. aborting...\n");
        abort();
    }

    
    int ipt_reads = 0;  //total reads that were analyzed
    int pass_reads = 0; //number of reads after fastp filtering
    
    char smpl_nm[MAX_LINE] = {0};   //name of samples
    get_smpl_nm(prnt_dir, smpl_nm); //obtain sample name from parent directory input
    
    if (!test_SM2_data) { //only read fastp output if not assessing test SM2 data
        get_fastp_out(prnt_dir, &ipt_reads, &pass_reads); //get filtering stats from fastp json output
    }
    get_profiles(prnt_dir, &mtrx, &algn[0], test_SM2_data); //store reactivity/eff depth of each transcript len in cotrans matrix
    
    
    if (mode == MULTI) {
        check_contiguity(&mtrx); //check that transcript lengths are contiguous from min len to max len
        set_enriched_lengths(&mtrx, incld_up2, excld_trmnl);
        calc_reads_algnd(&algn[0], &mtrx);
        
        //print csv files
        print_csv(prnt_dir, smpl_nm, "_reactivity.csv", &mtrx.vals, mtrx.tl[MIN], mtrx.tl[MAX], mtrx.nt[MIN], mtrx.nt[MAX]); //make reactivity csv
        print_csv(prnt_dir, smpl_nm, "_untreated.csv", &mtrx.unt, mtrx.tl[MIN], mtrx.tl[MAX], mtrx.nt[MIN], mtrx.nt[MAX]); //make mod ef depth csv
        print_csv(prnt_dir, smpl_nm, "_modified.csv", &mtrx.mod, mtrx.tl[MIN], mtrx.tl[MAX], mtrx.nt[MIN], mtrx.nt[MAX]); //make mod ef depth csv
        
        //print matrix columns
        columnize_mtrx(prnt_dir, smpl_nm, &mtrx);
        
        //print alignment details
        if (!test_SM2_data) {
            print_alignment_stats(prnt_dir, smpl_nm, &algn[0], &mtrx);
        }
        
    }
    
    mk_linebars(prnt_dir, smpl_nm, &mtrx); //format data for linebar plots
    
    if (rdat_config_provided) {
        mk_rdat(fp_config, &mtrx, mode); //generate rdat file
    }
}
















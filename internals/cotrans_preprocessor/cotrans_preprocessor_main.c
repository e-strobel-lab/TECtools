//
//  CTRL0_cotrans_preprocessor_main.c
//  
//
//  Created by Eric Strobel on 8/25/21.
//


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <sys/stat.h>
#include <math.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "cotrans_preprocessor_defs.h"
#include "cotrans_preprocessor_structs.h"
#include "../utils/io_management.h"
#include "../seq_utils/parse_fasta.h"
#include "../seq_utils/mk_fasta.h"
#include "./MLT_trgt_gen/mk_MLT_trgts.h"
#include "./SGL_trgt_gen/mk_SGL_target.h"
#include "./prcs_rds/MLT/prcs_MLT_cotrans.h"
#include "./run_script_gen/MLT/mk_MLT_run_script.h"
#include "../utils/debug.h"

extern int debug;			//flag to run debug mode
extern int debug_S2B_hash;	//seq2bin_hash-specific debug flag
extern struct testdata_3pEnd_vars testdata_3pEnd; //structure containing test data analysis variables

/* set_run_mode: set run mode using mode argument string */
int set_run_mode(char * mode_arg, int * run_mode, fastp_params * fastp_prms);

//functions to check that required input files were supplied for each mode
int check_make_fasta(int fa_nm_provided, int fa_sq_provided);
int check_make_3pEnd_targets(int fa_provided);
int check_standard_cotrans(int fq1_provided, int fq2_provided, int ends_provided, int fa_ref_provided, fastp_params fastp_prms);
int check_make_SM2_run_script(int config_provided);

//other functions
int get_target_len(FILE * ifp);


int main(int argc, char *argv[])
{
    FILE *fp_fasta = NULL;			//pointer to fasta file
    FILE *fp_3pEnd = NULL;			//pointer for 3' ends file
    FILE *fp_faRef = NULL;          //pointer for fasta reference file, used for TECprobe-SL
    FILE *fp_config_MLT = NULL;		//pointer to config file
    
    int run_mode = -1;				//run mode setting
    
    //fasta generation options
    char fasta_nm[MAX_LINE] = {0};
    char fasta_sq[MAX_LINE] = {0};
    
    //3' end target generation options
    int usr_end_len = 0;			//user-specified end target length; 0 if default
    int usr_min_len = 0;			//user-specified min transcript intermediate length, 0 if default
    int make_test_data = 0;			//flag to generate test data
    int randomize_end = 0;			//flag to randomized test data 3' ends
    int mltplir = 1;				//multiplier to make more test data reads;
    
    //variables to track option inputs
    int fa_nm_provided = 0;         //tracks number of names provided for fasta gen
    int fa_sq_provided = 0;         //tracks number of sequences provided for fasta gen
    int fq1_provided = 0;			//tracks number of read1 files provided
    int fq2_provided = 0;			//tracks number of read2 files provided
    int ends_provided = 0;			//tracks number of provided 3' end files
    int fa_ref_provided = 0;        //tracks number of provided fasta reference files
    int fastp_path_provided = 0;	//tracks number of fastp paths provided
    int fa_provided = 0;			//tracks number of fasta files provided
    int config_MLT_provided = 0;	//tracks number of config files provided
    
    int fgen_option_provided = 0;   //flag that fastq generation option was provided
    int ends_option_provided = 0;   //flag that 3' ends generation option was provided
    int prcs_option_provided = 0;   //flag that read processing option was provided
    int cnfg_option_provided = 0;   //flag that run script config option was provided
    
    names nm = {{{0}}};				//file and sample names
    fastp_params fastp_prms = {"fastp", -1, 0}; //parameters for fastp processing
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            //general arguments
            {"mode",         required_argument,  0,  'm'},  //set flag for run mode
            
            //fasta generation arguments
            {"seq-name",  required_argument, 0, 'n'},
            {"sequence",  required_argument, 0, 's'},
            
            //3' end target generation arguments
            {"fasta",				required_argument,	0,  'A'},	//fasta file input
            {"end-length",			required_argument,  0,  'E'},	//user-specified end target length
            {"min-length",			required_argument,  0,  'S'},	//minimum intermediate target length
            {"test-data",			no_argument,        0,  'T'},	//make test data
            {"multiplier",			required_argument,  0,  'U'},	//test data read multiplier
            {"end-randomization",	no_argument,  		0,  'R'},	//randomize test data 3' end
            //TODO: have the targets generation process make a hash table to assess unacceptable collisions ahead of analysis
            
            
            //read processing arguments
            {"read1",        required_argument,  0,  'i'},	//read 1 input
            {"read2",        required_argument,  0,  'I'},	//read 2 input
            {"3pEnd", 		 required_argument,  0,  'e'},	//3' end file input
            {"fasta-ref",    required_argument,  0,  'a'},  //fasta for single length data processing
            {"fastp-path",   required_argument,  0,  'p'},  //provide path to fastp executable
            {"debug",        no_argument,        0,  'd'},  //turn on debug mode
            {"testdata",     no_argument,        0,  't'},  //turn on test data read analysis mode
            {"limit",        no_argument,        0,  'l'},  //set limit on the number of reads to process
            
            //run script generation arguments
            {"config",       required_argument,  0,  'c'},	//config file input
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "m:n:s:A:E:S:TU:Ri:I:e:a:p:dtl:c:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
              
                
            /* set run mode and, if necessary, fastp processing mode */
            case 'm':
                set_run_mode(argv[optind-1], &run_mode, &fastp_prms);
                break;
                
                
            /* set sequence name for fasta generation */
            case 'n':
                fa_nm_provided++;
                strcpy(fasta_nm, argv[optind-1]);       //store fasta generation name
                fgen_option_provided = 1;
                break;
                
                
            /* set sequence for fasta generation */
            case 's':
                fa_sq_provided++;
                strcpy(fasta_sq, argv[optind-1]);       //store fasta generation sequence
                fgen_option_provided = 1;
                break;
              
                
            /* get fasta input file */
            case 'A':
                fa_provided++;                         //count fasta files provided
                get_file(&(fp_fasta), argv[optind-1]); //set file pointer to input fasta file
                ends_option_provided = 1;
                break;
                
                
                
            /* get user-specified 3' end target length */
            case 'E':
                usr_end_len = atoi(argv[optind-1]);	   //set 3' end target length
                if (usr_end_len > MAX_END_LEN) {	   //check that usr_end_len is less than the allowed max
                    printf("cotrans_preprocessor_main: error - end target length (%d) is greater than allowed maximum (%d). aborting...\n", usr_end_len, MAX_END_LEN);
                    abort();
                }
                ends_option_provided = 1;
                break;
                
                
            /* set user-specified minimum transcript intermediate length */
            case 'S':
                usr_min_len = atoi(argv[optind-1]); //set minimum transcript intermediate length
                
                //check that usr_min_len is above the allowed minimum
                //the comparison is to DFLT_MIN+1 because the minimum transcript length will
                //eventually be set as usr_min_len-1 in the set_min_len function (mk_MLT_trgts.c)
                if (usr_min_len < DFLT_MIN+1) {
                    printf("cotrans_preprocessor_main: error - shortest transcript length(%d) is less than allowed minimum (%d). aborting...\n", usr_min_len, DFLT_MIN+1);
                    abort();
                }
                ends_option_provided = 1;
                break;
                
                
            /* turn test data generation on */
            case 'T':
                make_test_data = 1;
                ends_option_provided = 1;
                break;
            
                
            /* set read multiplier value */
            case 'U':
                mltplir = atoi(argv[optind-1]);
                ends_option_provided = 1;
                break;
            
                
            /* turn test data 3' end randomization on */
            case 'R':
                randomize_end = 1;
                ends_option_provided = 1;
                break;
            
                
            /* get read 1 input filename */
            case 'i':
                fq1_provided++;							//count read 1 fastq files provided
                strcpy(nm.file[READ1], argv[optind-1]);	//store read 1 filename
                prcs_option_provided = 1;
                break;
            
                
            /* get read 2 input filename */
            case 'I':
                fq2_provided++;							//count read 2 fastq files provided
                strcpy(nm.file[READ2], argv[optind-1]);	//store read 2 filename
                prcs_option_provided = 1;
                break;
            
            
                
            /* get 3' end targets input file */
            case 'e':
                ends_provided++;						//count 3' end targets files provided
                strcpy(nm.ends, argv[optind-1]);		//store 3' end targets filename
                get_file(&(fp_3pEnd), argv[optind-1]);	//set file pointer to input 3' ends target file
                prcs_option_provided = 1;
                break;
            
            case 'a':
                fa_ref_provided++;                      //count fasta reference files provided
                get_file(&(fp_faRef), argv[optind-1]);  //set file pointer to input fasta reference target file
                prcs_option_provided = 1;
                break;
                
            /* set path to fastp executable */
            case 'p':
                fastp_path_provided++;
                strcpy(fastp_prms.path, argv[optind-1]);
                prcs_option_provided = 1;
                break;
            
                
            /* turn on debug mode */
            case 'd':
                debug = 1;			//flag for general debug mode
                debug_S2B_hash = 1; //flag for seq2bin hash debug mode //TODO: make separate option?
                break;
            
                
            /* turn on test data read analysis mode */
            case 't':
                testdata_3pEnd.run = 1;
                prcs_option_provided = 1;
                break;
            
                
            /* set read processing limit */
            case 'l':
                fastp_prms.limit = atoi(argv[optind-1]);
                prcs_option_provided = 1;
                break;
            
                
            /* get config file */
            case 'c':
                config_MLT_provided++;						//count config files provided
                get_file(&(fp_config_MLT), argv[optind-1]); //set file pointer to input config file
                cnfg_option_provided = 1;
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
    
    
    /*********** branch to user-specified run mode  ***********/
    switch (run_mode) {

        /****************** fasta generation ******************/
        case MK_FASTA:
            if (ends_option_provided || prcs_option_provided || cnfg_option_provided) {
                printf("cotrans_preprocessor_main: error - run mode set to MK_FASTA but at least one argument unrelated to fasta generation was provided. aborting...");
                abort();
            }
            
            check_make_fasta(fa_nm_provided, fa_sq_provided);
            
            mk_fasta_file(fasta_nm, fasta_sq);
            
            break;
        /*************** end of fasta generation **************/
            
        
            
        /************** 3' end target generation **************/
        case MK_3pEND_TARGETS:
            
            //check that unnecessary arguments were not provided
            if (fgen_option_provided || prcs_option_provided || cnfg_option_provided) {
                printf("cotrans_preprocessor_main: error - run mode set to MK_3pEND_TARGETS but at least one argument unrelated to multi-length target generation was provided. aborting...");
                abort();
            }
            
            //check that required arguments were provided correctly
            check_make_3pEnd_targets(fa_provided);
            
            //start 3' end target generation
            mk_MLT_trgts(fp_fasta, usr_min_len, usr_end_len, make_test_data, mltplir, randomize_end);
            
            break;
        /*********** end of 3' end target generation **********/
        
        
        
        /************** single target generation **************/
        case MK_SINGLE_TARGET:
            //check that unnecessary arguments were not provided
            if (fgen_option_provided || prcs_option_provided || cnfg_option_provided) {
                printf("cotrans_preprocessor_main: error - run mode set to MK_SINGLE_TARGET but at least one argument unrelated to multi-length target generation was provided. aborting...");
                abort();
            }
            
            //check that required arguments were provided correctly
            //NOTE: the minimum required arguments are the same for
            //multi-length and single targets, so it is currently valid
            //to use check_make_3pEnd_targets to assess whether required
            //arguments were provided when making single targets.
            check_make_3pEnd_targets(fa_provided);
            
            mk_SGL_target(fp_fasta);
            
        break;
        /*********** end of single target generation **********/
            
        /************ sequencing read preprocessing ***********/
        case PRCS_READS:
            
            //check that unnecessary arguments were not provided
            if (fgen_option_provided || ends_option_provided || cnfg_option_provided) {
                printf("cotrans_preprocessor_main: error - run mode set to PRCS_READS but at least one argument unrelated to sequencing read processing was provided. aborting...");
                abort();
            }
            
            //standard cotrans read processing
            if (fastp_prms.mode == MULTI || fastp_prms.mode == SINGLE) {
                
                //check that required arguments were provided correctly
                check_standard_cotrans(fq1_provided, fq2_provided, ends_provided, fa_ref_provided, fastp_prms);
                if (testdata_3pEnd.run) {
                    check_testdata_ipt(&nm);
                }
                
                if (fastp_prms.mode == SINGLE) {       //when processing TECpobe-SL data
                    //get length of target from fasta reference
                    //note that the only purpose of the fasta reference is to determine the
                    //length of the target sequence, which is needed for run script generation
                    nm.len = get_target_len(fp_faRef);
                }
                
                //start sequencing read processing
                prcs_MLT_cotrans(&nm, fp_3pEnd, fastp_prms);
            }
            
            break;
        /******** end of sequencing read preprocessing ********/
        
            
        /**************** run script generation ***************/
        case MK_RUN_SCRIPT:
            
            //check that unnecessary arguments were not provided
            if (fgen_option_provided || ends_option_provided || prcs_option_provided) {
                printf("cotrans_preprocessor_main: error - run mode set to MK_RUN_SCRIPT but at least one argument unrelated to run script generation was provided. aborting...");
                abort();
            }
            
            //check that required arguments were provided correctly
            check_make_SM2_run_script(config_MLT_provided);
            
            //start run script generation
            mk_MLT_run_script(fp_config_MLT);
            
            break;
        /************* end of run script generation ***********/
            
        //no run mode specified
        case -1:
            printf("cotrans_preprocessor_main: run mode was not specified. aborting...\n");
            abort();
            break;
            
        //unexpected run mode
        default:
            printf("cotrans_preprocessor_main: unexpected run mode. aborting...\n");
            abort();
            break;
    }
}

/* set_run_mode: set run mode using mode argument string */
int set_run_mode(char * mode_arg, int * run_mode, fastp_params * fastp_prms)
{
    if (!strcmp(mode_arg, "MAKE_FASTA")) {
        *run_mode = MK_FASTA;
        return 1;
    } else if (!strcmp(mode_arg, "MAKE_3pEND_TARGETS")) {
        *run_mode = MK_3pEND_TARGETS;
        return 1;
    } else if (!strcmp(mode_arg, "MAKE_SINGLE_TARGET")) {
        *run_mode = MK_SINGLE_TARGET;
        return 1;
    } else if (!strcmp(mode_arg, "PROCESS_MULTI")) {
        *run_mode = PRCS_READS;
        fastp_prms->mode = MULTI;
        return 1;
    } else if (!strcmp(mode_arg, "PROCESS_SINGLE")) {
        *run_mode = PRCS_READS;
        fastp_prms->mode = SINGLE;
        return 1;
    } else if (!strcmp(mode_arg, "MAKE_RUN_SCRIPT")) {
        *run_mode = MK_RUN_SCRIPT;
        return 1;
    } else {
        printf("set_run_mode: error - unrecognized mode. please use MAKE_FASTA, MAKE_3pEND_TARGETS, MAKE_SINGLE_TARGET, PROCESS_MULTI, or PROCESS_SINGLE to indicate run mode.\n");
        abort();
    }
}

/* check_make_fasta: check that required inputs for fasta generation were provided */
int check_make_fasta(int fa_nm_provided, int fa_sq_provided)
{
    if (fa_nm_provided != 1) {
        printf("cotrans_preprocessor_main: error incorrect number (%d) of fasta output names provided\n", fa_nm_provided);
        abort();
    }
    
    if (fa_sq_provided != 1) {
        printf("cotrans_preprocessor_main: error incorrect number (%d) of fasta output sequences provided\n", fa_sq_provided);
        abort();
    }
    
    return 1;
}

/* check_make_3pEnd_targets: check that required inputs for 3' end target gen were provided */
int check_make_3pEnd_targets(int fa_provided)
{
    //check that only one fasta file was provided
    if (fa_provided != 1) {
        printf("cotrans_preprocessor_main: error - incorrect number (%d) of fasta files provided\n", fa_provided);
        abort();
    }
    
    return 1;
}

/* check_standard_cotrans: check that required inputs for standarrd cotrans preprocessing were provided */
int check_standard_cotrans(int fq1_provided, int fq2_provided, int ends_provided, int fa_ref_provided, fastp_params fastp_prms)
{
    //check that only one read1 fastq file was provided
    if (fq1_provided != 1) {
        printf("cotrans_preprocessor_main: error - incorrect number (%d) of read 1 fastq files provided\n", fq1_provided);
        abort();
    }
    
    //check that only one read2 fastq file was provided
    if (fq2_provided != 1) {
        printf("cotrans_preprocessor_main: error - incorrect number (%d) of read 2 fastq files provided\n", fq2_provided);
        abort();
    }
    
    //check that only one 3' ends target file was provided if experiment is multi-length
    if (ends_provided != 1 && fastp_prms.mode == MULTI) {
        printf("cotrans_preprocessor_main: error - incorrect number (%d) of 3' end targets files provided\n", ends_provided);
        abort();
    }
    
    //throw error if 3' ends target file was provided but multi-length option was not set
    if (ends_provided && fastp_prms.mode != MULTI) {
        printf("cotrans_preprocessor_main: error - 3' ends targets provided but run mode was not set to MULTI\n");
        abort();
    }
    
    //check that only one fasta reference file was provided if experiment is single length
    if (fa_ref_provided != 1 && fastp_prms.mode == SINGLE) {
        printf("cotrans_preprocessor_main: error - incorrect number (%d) of fasta reference files provided\n", fa_ref_provided);
        abort();
    }
    
    //throw error if fasta reference file was provided but single length option was not set
    if (fa_ref_provided && fastp_prms.mode != SINGLE) {
        printf("cotrans_preprocessor_main: error - fasta reference provided but run mode was not set to SINGLE\n");
        abort();
    }
    
    //check that processing mode was set
    if (fastp_prms.mode == -1) {
        printf("cotrans_preprocessor_main: error - fastp processing mode was not set\n");
        abort();
    }

    return 1;
}

/* check_make_SM2_run_script: check that required inputs for shapemapper 2 run script gen were provided */
int check_make_SM2_run_script(int config_MLT_provided)
{
    //check that only one config file was provided
    if (config_MLT_provided != 1) {
        printf("cotrans_preprocessor_main: error - incorrect number (%d) of config files provided. Aborting program...\n", config_MLT_provided);
        abort();
    }
        
    return 1;
}

/* get_target_len: get length of aligment target for single length data processing.
 target length is used when naming output files and in the config file */
int get_target_len(FILE * ifp)
{
    char line[MAX_LINE] = {0}; //array to store fasta file lines
    int len = 0;               //length of target
    
    static const char ldr[21]  = "atggccttcgggccaa"; //leader sequence, masked by lowercase
    
    get_line(&line[0], ifp);   //get first line of fasta file
    if (line[0] != '>') {
        printf("get_target_len: error - unexpected format for fasta file. aborting...");
        abort();
    }
    
    get_line(&line[0], ifp);   //get second line of fasta file
    len = strlen(line);        //get length of target sequence
    
    if (strstr(line, ldr) != NULL) { //check if target contains leader sequence
        len -= strlen(ldr);          //subtract length of leader from target length
    }
    
    return len;
}


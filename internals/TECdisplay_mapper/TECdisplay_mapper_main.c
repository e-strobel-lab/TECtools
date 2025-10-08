//
//  TECDisplay_mapper_main.c
//  
//
//  Created by Eric Strobel on 6/21/22.
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

#include "./TECdisplay_mapper_structs.h"
#include "./TECdisplay_mapper_defs.h"

#include "../utils/io_management.h"
#include "../utils/debug.h"

#include "../variant_maker/vmt_suffix.h"

#include "./map_reads/map_reads.h"

extern int debug;            //flag to run debug mode
extern int debug_S2B_hash;   //seq2bin_hash-specific debug flag

int set_run_mode(char * mode_arg, int * run_mode);
int set_min_qscore(char * min_Qscore, char * val2set);
int check_options(int fq1_provided, int fq2_provided, int trgs_provided);
int check_testdata_options(int fq1_provided, int fq2_provided, int trgs_provided, fastp_params * fastp_prms);

int main(int argc, char *argv[])
{
    extern char vmt_suffix[4]; //variant maker targets file suffix
    
    testdata_vars testdata = {0};   //structure for tracking testdata metrics
    FILE *fp_trgs = NULL;           //pointer to targets file
    int run_mode = -1;              //run mode setting

    //variables to track option inputs
    int fq1_provided = 0;           //tracks number of read1 files provided
    int fq2_provided = 0;           //tracks number of read2 files provided
    int trgs_provided = 0;          //tracks number of targets files provided
    int fastp_path_provided = 0;    //tracks number of fastp paths provided
    int min_vbase_qscore_set = 0;   //tracks if min variable base qscore was set
    int min_cbase_qscore_set = 0;   //tracks if min constant base qscore was set
    
    TDSPLY_names nm = {{{0}}};                  //file and sample names
    fastp_params fastp_prms = {"fastp", -1, 0}; //parameters for fastp processing
    
    char min_qscore[2] = {0};     //array of quality score thresholds
    min_qscore[Q_VARIABLE] = '5'; //initialize variable base qscore threshold to 20 (ascii: 5);
    min_qscore[Q_CONSTANT] = '!'; //initialize constant base qscore threshold to  0 (ascii: !);
    
    char * file_suffix = {0};        //pointer to file suffix
    int trgt_ftype = FILE_TYPE_INIT; //target file type
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"mode",          required_argument,  0,  'm'},  //mode
            {"read1",         required_argument,  0,  'i'},  //read 1 input
            {"read2",         required_argument,  0,  'I'},  //read 2 input
            {"targets",       required_argument,  0,  't'},  //targets file input
            {"out-name",      required_argument,  0,  'o'},  //output file name
            {"fastp-path",    required_argument,  0,  'p'},  //provide path to fastp executable
            {"qual-variable", required_argument,  0,  'v'},  //min quality for variable bases
            {"qual-constant", required_argument,  0,  'c'},  //min quality for constant bases
            {"debug",         no_argument,        0,  'd'},  //turn on debug mode
            {"limit",         required_argument,  0,  'l'},  //set limit on the number of reads to process
            
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "m:i:I:t:o:p:v:c:dl:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
            /* set run mode and, if necessary, fastp processing mode */
            case 'm':
                set_run_mode(argv[optind-1], &run_mode);
                break;
                
            /* get read 1 input filename */
            case 'i':
                fq1_provided++;                            //count read 1 fastq files provided
                strcpy(nm.file[READ1], argv[optind-1]);    //store read 1 filename
                break;
            
                
            /* get read 2 input filename */
            case 'I':
                fq2_provided++;                            //count read 2 fastq files provided
                strcpy(nm.file[READ2], argv[optind-1]);    //store read 2 filename
                break;
            
                
            /* get targets input file */
            case 't':
                trgs_provided++;                                        //count targets files provided
                strcpy(nm.trgs, argv[optind-1]);                        //store targets filename
                file_suffix = get_sample_name(nm.trgs, nm.trgs_prefix); //get targets sample name
                
                if (!strcmp(file_suffix, vmt_suffix)) { //check and set barcoded targets file type
                    trgt_ftype = VMT_FILE;
                } else if (!strcmp(file_suffix, "fasta") || !strcmp(file_suffix, "fa")) {
                    trgt_ftype = FASTA_FILE;
                } else {
                    printf("TECdisplay_mapper: error - unrecognized targets file type (.%s). aborting...\n", file_suffix);
                    abort();
                }
                
                get_file(&(fp_trgs), argv[optind-1]);       //set file pointer to input target file
                break;
                
            /* get output filename */
            case 'o':
                strcpy(nm.out_nm, argv[optind-1]);    //store output filename
                break;
                
            /* set path to fastp executable */
            case 'p':
                fastp_path_provided++;
                strcpy(fastp_prms.path, argv[optind-1]);
                break;
                
            /* set variable base qscore threshold */
            case 'v':
                if (!min_vbase_qscore_set) {
                    set_min_qscore(&min_qscore[Q_VARIABLE], argv[optind-1]);
                } else {
                    printf("main: error - more than one minimum variable base quality score was supplied. aborting...\n");
                    abort();
                }
                break;
                
            /* set constant base qscore threshold */
            case 'c':
                if (!min_cbase_qscore_set) {
                    set_min_qscore(&min_qscore[Q_CONSTANT], argv[optind-1]);
                } else {
                    printf("main: error - more than one minimum constant base quality score was supplied. aborting...\n");
                    abort();
                }
                break;
                
            /* turn on debug mode */
            case 'd':
                debug = 1;          //flag for general debug mode
                debug_S2B_hash = 1; //flag for seq2bin hash debug mode
                break;
                            
                
            /* set read processing limit */
            case 'l':
                fastp_prms.limit = atoi(argv[optind-1]);
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
    
    printf("min qscore for variable bases: %2d (%c)\n", min_qscore[Q_VARIABLE]-'!', min_qscore[Q_VARIABLE]);
    printf("min qscore for constant bases: %2d (%c)\n", min_qscore[Q_CONSTANT]-'!', min_qscore[Q_CONSTANT]);
        
    if (run_mode == MAP_SEQ_READS) {                                   //run MAP_SEQ_READS mode
        check_options(fq1_provided, fq2_provided, trgs_provided);      //check that correct options were supplied
        map_reads(&nm, fp_trgs, trgt_ftype, min_qscore, fastp_prms, &testdata, MAP_SEQ_READS); //map sequencing reads
        
                
    } else if (run_mode == MAP_TEST_DATA) {                               //run MAP_TEST_DATA mode
        check_testdata_options(fq1_provided, fq2_provided, trgs_provided, &fastp_prms); //check options
        testdata.run = 1;                                                 //turn testdata run flag on
        map_reads(&nm, fp_trgs, trgt_ftype, min_qscore, fastp_prms, &testdata, MAP_TEST_DATA   ); //map test data reads to target
        
    } else {
        printf("TECDisplay_mapper_main: error - unrecognized run mode. aborting...\n");
        abort();
    }
}

/* set_run_mode: set run mode using mode argument string */
int set_run_mode(char * mode_arg, int * run_mode)
{
    if (!strcmp(mode_arg, "MAP_TEST_DATA")) {
        *run_mode = MAP_TEST_DATA;
        return 1;
    } else if (!strcmp(mode_arg, "MAP_SEQ_READS")) {
        *run_mode = MAP_SEQ_READS;
        return 1;
    } else {
        printf("set_run_mode: error - unrecognized mode. please use MAP_SEQ_READS or MAP_TEST_DATA to indicate run mode.\n");
        abort();
    }
}

/* set_min_qscore: check that minimum qscore input is valid and set min_qscore variable */
int set_min_qscore(char * min_qscore, char * val2set)
{
    int i = 0;                 //general purpose index
    int len = strlen(val2set); //length of the input string
    int val = 0;               //value of the input string
    
    //check that val2set string is 1 or 2 characters long
    if (len < 1 || len > 2) {
        printf("set_min_qscore: error - minimum qscore input should be a value >=0 and <=41 but has a string length of %d. aborting...\n", len);
        abort();
    }
    
    //check that val2set is entirely composed of digits
    for (i = 0; i < len; i++) {
        if (!isdigit(val2set[i])) {
            printf("set_min_qscore: error - minimum qscore input should be a value >=0 and <=41 but is not solely composed of digits. aborting...\n");
            abort();
        }
    }
    
    val = atoi(val2set);       //convert input string to integer
    if (val < 0 || val > 41) { //if minimum qscore is out of the permissible range
        printf("set_min_qscore: error - minimum qscore values must be >=0 and <=41\n");
        abort();
    } else {
        *min_qscore = '!' + (char)val;
    }
    
    return 1;
}

/* check_options: check that required inputs for sequencing read mapping were provided */
int check_options(int fq1_provided, int fq2_provided, int trgs_provided)
{
    //check that only one read1 fastq file was provided
    if (fq1_provided != 1) {
        printf("TECdisplay_mapper_main: error - incorrect number (%d) of read 1 fastq files provided\n", fq1_provided);
        abort();
    }
    
    //check that only one read2 fastq file was provided
    if (fq2_provided != 1) {
        printf("TECdisplay_mapper_main: error - incorrect number (%d) of read 2 fastq files provided\n", fq2_provided);
        abort();
    }
    
    //check that only one target files was provided
    if (trgs_provided != 1) {
        printf("TECdisplay_mapper_main: error - incorrect number (%d) of targets files provided\n", trgs_provided);
        abort();
    }

    return 1;
}

/* check_testdata_options: check that required inputs for test data analysis were provided */
int check_testdata_options(int fq1_provided, int fq2_provided, int trgs_provided, fastp_params * fastp_prms)
{
    //check that no fastq files were provided
    if (fq1_provided || fq2_provided) {
        printf("TECdisplay_mapper_main: error - one or more fastq files were provided when running MAP_TEST_DATA mode. aborting...\n");
        abort();
    }
    
    //check that only one target files was provided
    if (trgs_provided != 1) {
        printf("TECdisplay_mapper_main: error - incorrect number (%d) of targets files provided\n", trgs_provided);
        abort();
    }
    
    //check that read processing limit was not set
    if (fastp_prms->limit) {
        printf("TECdisplay_mapper_main: error - read processing limit was set when running MAP_TEST_DATA mode. aborting...");
        abort();
    }

    return 1;
}

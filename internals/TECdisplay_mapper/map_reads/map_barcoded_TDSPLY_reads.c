//
//  map_barcoded_TDSPLY_reads.c
//  
//
//  Created by Eric Strobel on 10/14/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../cotrans_preprocessor/MUX_trgt_gen/mk_MUX_trgts.h"
#include "../../cotrans_preprocessor/prcs_rds/MUX/map_barcoded_targets.h"
#include "../../cotrans_preprocessor/prcs_rds/MUX/get_brcd_str.h"
#include "../../cotrans_preprocessor/prcs_rds/MUX/testdataMUX_analysis.h"

#include "../TECdisplay_mapper_defs.h"
#include "../TECdisplay_mapper_structs.h"

#include "../../utils/io_management.h"
#include "../../utils/debug.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "./UNV/call_fastp_TDSPLY.h"
#include "./UNV/prcs_chnl_TDSPLY.h"

#include "../../seq_utils/revcomp.h"
#include "../../seq_utils/basemap.h"
#include "../../seq_utils/mapping_metrics.h"

#include "./map_expected/mk_output_files.h"
#include "./map_expected/print_navigator_template.h"

#include "../testdata_analysis/mk_BC_TDSPLY_test_data.h"
//#include "../testdata_analysis/assess_TDSPLY_test_data.h"

#include "map_barcoded_TDSPLY_reads.h"

extern int debug;

/* prcs_barcoded_TDSPLY_reads: coordinates targets parsing, fastp processing, read mapping, and output file generation */
int prcs_barcoded_TDSPLY_reads(TDSPLY_names * nm, FILE * fp_trgs, int trgt_ftype, char * minQ, fastp_params fastp_prms, testdata_vars * testdata, int mode)
{
    //TODO: need system to ignore linker (if present) during mapping verification.
    //TODO: check how variant_maker handles upper and lower case letters when generating targets
    
    extern testdata_MUX_vars testdata_MUX; //variables for testdata analysis - using MUX code here
    testdata_MUX.run = testdata->run;      //hand off run variable from testdata to testdata_MUX
    
    FILE *ifp = NULL;           //pointer for merged fastq file
    mapping_metrics met = {0};  //read processing metrics storage
    
    init_chnl_mtrcs_mem(&met, TDSPLY_CHANNEL_MAX); //initialize channel tracking memory
    
    /************* targets initialization and parsing **************/
    target *refs = NULL;          //pointer for array of reference targets
    opt_ref *ref_val = NULL;      //pointer for array of optional reference target structures
    
    compact_target *ctrg = NULL;  //pointer for array of hash table targets
    opt_BC *BC_val = NULL;        //pointer for array of optional target structures
    
    target_params trg_prms = {0}; //structure for storing target parameters
    TDSPLY_fasta wt = {0};        //storage for wt sequence information
    
    int ctrg_cnt = 0;             //number of compact targets stored
    int clcd_ctrg_cnt = 0;        //calculated number of barcode targets
    
    char line[MAX_LINE+1] = {0};  //array to store line
    
    int ret = 0; //variable for storing snprintf return value
    
    //TODO: 99% sure that this can be deleted
    //allocate memory for reference targets
    if ((refs = calloc(MAXREF, sizeof(*refs))) == NULL) {
        printf("map_reads: error - reference target memory allocation failed\n");
        return 1;
    }
    
    if ((ref_val = calloc(MAXREF, sizeof(*ref_val))) == NULL) {
        printf("map_reads: error - reference target value memory allocation failed\n");
        return 1;
    }
    
    //determine expected target count
    if (trgt_ftype == FASTA_FILE) {
        
        //iterate through file to count number of targets
        while (get_line(line, fp_trgs)) { //until all lines have been read
            if (line[0] == '>') {         //if reading first line of fasta entry
                get_line(line, fp_trgs);  //get the second line of the fasta entry
                trg_prms.xpctd++;         //increment expected barcode count
            } else {
                printf("prcs_MUX_cotrans: error - unexpected format for targets fasta file. aborting...\n");
                abort();
            }
        }
        
        fclose(fp_trgs);                //close targets file
        get_file(&(fp_trgs), nm->trgs); //re-open targets file
        
    } else {
        printf("map_barcoded_TDSPLY_reads: expected targets as fasta file. aborting...\n");
        abort();
    }
    
    clcd_ctrg_cnt = trg_prms.xpctd * 129; //caculate expected number of BC targs, including single subs and indels
    
    //allocate memory for compact targets
    if ((ctrg = calloc(clcd_ctrg_cnt, sizeof(*ctrg))) == NULL) {
        printf("map_barcoded_TDSPLY_reads: error - target memory allocation failed\n");
        return 1;
    }
    
    //allocate memory for barcode target optional values
    if ((BC_val = calloc(clcd_ctrg_cnt, sizeof(*BC_val))) == NULL) {
        printf("map_barcoded_TDSPLY_reads: error - target value memory allocation failed\n");
        return 1;
    }
    
    printf("\nProcessing targets file that contains %d targets\n\n", trg_prms.xpctd);
    
    //generate barcode targets
    ctrg_cnt = mk_MUX_trgts(refs, ref_val, ctrg, BC_val, fp_trgs, trgt_ftype, &trg_prms, clcd_ctrg_cnt, &wt, TDSPLY);
    met.srcTrgs = trg_prms.t_cnt; //record number of source barcodes
    met.targets = ctrg_cnt;       //record number of targets generated
    
    /********** end of targets initialization and parsing **********/

    /********* hash table initialization and construction **********/
    compact_h_node **htbl = NULL;                  //hash table root
    compact_h_node_bank hn_bank = {NULL, NULL, 0}; //bank for hash table nodes
    
    hn_bank.count = 0; //initialize hash table node count to zero
    
    //allocate TABLE_SIZE hash table node pointers
    if ((htbl = calloc(TABLE_SIZE, sizeof(*htbl))) == NULL) {
        printf("map_barcoded_TDSPLY_reads: error - hash table memory allocation failed. aborting...\n");
        abort();
    }
    
    //allocate BLOCK_SIZE hash table nodes
    if ((hn_bank.chn = calloc(BLOCK_SIZE, sizeof(*(hn_bank.chn)))) == NULL) {
        printf("map_barcoded_TDSPLY_reads: error - hash table node memory allocation failed. aborting...\n");
        abort();
    }
    
    compact_h_node_bank *crrnt_hn_bank = &hn_bank;                     //pointer for handling hash node bank
    mk_htbl_MUX(htbl, crrnt_hn_bank, ctrg, ctrg_cnt, &trg_prms, &met); //generate barcode target hash table
    /****** end of hash table initialization and construction ******/

    
    /*************** testdata generation *****************/
     if (mode == MAP_TEST_DATA) {
        mk_BC_TDSPLY_test_data(nm, ctrg, &trg_prms, ctrg_cnt); //generate test data
     }    
    /*********** end of test data generation *************/
    
    
    /***************** obtain sample name prefix *******************/
    
    get_sample_name(nm->file[READ1], nm->smpl[READ1]); //get read 1 sample name
    get_sample_name(nm->file[READ2], nm->smpl[READ2]); //get read 2 sample name
    ret = snprintf(nm->mrg, MAX_LINE, "%s_merged", nm->smpl[READ1]); //construct merged read sample name
    if (ret >= MAX_LINE || ret < 0) {
        printf("map_barcoded_TDSPLY_reads: error - error when constructing output file name. aborting...\n");
        abort();
    }
    
    /************** end of obtain sample name prefix ***************/
    
    /************* process sequencing reads **************/
    
    mk_out_dir("processed");           //make directory for output files
    call_fastp_TDSPLY(nm, fastp_prms); //fastp pre-processing todo pass by reference?
    
    char processed_file[MAX_LINE+1]; //array for storing fastp-processed file name
    char gunzip[MAX_LINE+1];         //array for storing gunzip command
    
    //construct merged read file name
    ret = snprintf(processed_file, MAX_LINE, "./processed/%s.fq", nm->mrg);
    if (ret >= MAX_LINE || ret < 0) {
        printf("map_barcoded_TDSPLY_reads: error - error when constructing merged read file name. aborting...\n");
        abort();
    }
    
    //construct gunzip command
    ret = snprintf(gunzip, MAX_LINE, "gunzip %s", processed_file);
    if (ret >= MAX_LINE || ret < 0) {
        printf("map_barcoded_TDSPLY_reads: error - error when constructing gunzip command. aborting...\n");
        abort();
    }
    
    system(gunzip); //gunzip merged reads file
    
    //open merged read file
    if ((ifp = fopen(processed_file, "r")) == NULL) { //open merged reads file
        printf("map_barcoded_TDSPLY_reads: error - could not open %s as a merged read file. Aborting program...\n", processed_file);
        abort();
    }
    
    //map reads
    map_barcoded_TDSPLY_reads(ifp, htbl, ctrg, &trg_prms, &met, mode);
    
    //count non-redundant targets with >=1 mapped read
    trg_prms.mapped2 = count_matched_compact_targets(ctrg, &trg_prms, ctrg_cnt);
    
    //print output file
    print_output(ctrg, &trg_prms, nm, fastp_prms.mode);
     
    //close merged read file
    if ((fclose(ifp)) == EOF) {
        printf("map_barcoded_TDSPLY_reads: error - error occurred when closing merged read file. Aborting program...\n");
        abort();
    }
    
    /********** end of process sequencing reads **********/
    
    /************ print metrics and clean up *************/
    
    print_metrics(&trg_prms, &met, nm); //print mapping metrics
    
    if (testdata->run) { //print test data analysis
        print_MUX_testdata_analysis(&met, ctrg, TDSPLY);
    }
    
    //compress merged fastq file
    char gzip[MAX_LINE+1]; //array to store gzip command
    
    //construct gzip command
    ret = snprintf(gzip, MAX_LINE, "gzip %s", processed_file);
    if (ret >= MAX_LINE || ret < 0) {
        printf("map_reads: error - error when constructing gzip command. aborting...\n");
        abort();
    }
    system(gzip); //gzip merged read file
    
    /*********** end of print metrics clean up ***********/
    
    //free memory
    free(refs);
    free(ctrg);
    free(ref_val);
    free(BC_val);
    free(htbl);
    free(hn_bank.chn);
    
    return 1;
}

/* map_barcoded_TDSPLY_reads: map reads to user-supplied fasta barcoded targets */
void map_barcoded_TDSPLY_reads(FILE *ifp, compact_h_node **htbl, compact_target *ctrg, target_params * trg_prms, mapping_metrics * met, int mode)
{
    extern int debug;                      //flag to turn on debug mode
    extern testdata_MUX_vars testdata_MUX; //variables for testdata analysis - using MUX code here
    
    //anchor variables
    const char SC1_anchr[15] = "GGCCTTCGGGCCAA";   //structure cassette 1 sequence; anchor for finding RNA 5' end
    int anchr_len = strlen(SC1_anchr);             //length of SC1 sequence
    
    //general variables
    int i = 0;       //general purpose index
    int j = 0;       //general purpose index
    int proceed = 1; //flag to indicate that read processing should proceed
                
    //read variables
    char read[FQ_LINES][MAX_LINE+1] = {{0}};  //array for storing all four lines of fastq entry
    char read_rc[MAX_LINE+1] = {0};           //array for storing reverse complement of read
    
    //mapping variables
    char brcd_str[MAX_BARCODE_LEN+1] = {0};    //barcode string
    char rc_brcd_str[MAX_BARCODE_LEN+1] = {0}; //reverse complement of barcode string
    
    char * anchr_pt = NULL; //pointer to the location of the SC1 hairpin anchor sequence
    char * end5p = NULL;    //pointer to 5' end of the RNA target sequence (after SC1 hairpin)
    
    compact_target * crnt_ref_trg = NULL; //pointer to reference target of mapped target
    compact_target * crnt_mpd_trg = NULL; //pointer to mapped target
    
    opt_BC * BC_val = NULL; //pointer to barcode target optional values
    
    int ref_indx = 0; //index for reference targets
    int rd_indx = 0;  //index for end5p array
    int mtch = 0;     //flag that read is a match to target
    
    //channel variables
    int  channel = -1;       //stores channel barcode code
    int  chnl_mtch_typ = -1; //stores channel match type (full/partial)
    
    /* testdata variables */
    char id_line[MAX_LINE+1] = {0}; //array to store testdata id line (fastq line 1)
    char * td_trg_id = NULL;        //pointer to expected target id in test data id line
    int crnt_mut_cd = -1;           //current mutation code
    
    for (i = 0; proceed; i++) {
        
        /*** zero critical values before each loop iteration ***/
        for (j = 0; j < FQ_LINES; j++) { //initialize fastq arrays to 0
            read[j][0] = '\0';
        }
        
        read_rc[0] = '\0';    //initialize read revcomp array to 0
        
        anchr_pt = NULL;      //initialize anchor point pointer to NULL
        end5p = NULL;         //initialize 5' end pointer to NULL
        brcd_str[0] = '\0';   //initialize barcode string to 0
        
        crnt_ref_trg = NULL;  //initialize ref target pointer to NULL
        crnt_mpd_trg = NULL;  //initialize mapped target pointer to NULL
        
        channel = -1;         //initialize channel value to -1
        chnl_mtch_typ = -1;   //initialize channel match type to -1
        
        BC_val = NULL;        //initialize barcode vals pointer to NULL
        mtch = 0;             //initialize mtch to 0 (done below in loop as well; here for completeness)
        
        id_line[0] = '\0';    //initialize testdata id line to 0
        td_trg_id = NULL;     //initialize target id pointer to NULL
        crnt_mut_cd = -1;     //initialize current mutation code to -1
        
        /* copy fastq file lines to read array */
        for (j = 0; j < FQ_LINES && proceed; j++) {
            if (!get_line(&read[j][0], ifp)) { //get merged read fastq line, test success of get_line
                //test failed. this should only occur at end of file and get_line will
                //throw an error if an unexpected input line is encountered
                proceed = 0;
            }
        }
        
        if (proceed) {
            met->reads_processed++; //track total number of reads processed
            
            get_brcd_str(brcd_str, &read[LINE1][0]);            //get barcode from read1 id
            reverse_complement(rc_brcd_str, brcd_str, REVCOMP); //revcomp barcode string
                        
            crnt_ref_trg = map_brcd(brcd_str, rc_brcd_str, htbl, &crnt_mpd_trg, met); //map barcode using hash table
            
            if (crnt_ref_trg != NULL) { //if barcode mapped
                
                BC_val = (opt_BC *)(crnt_ref_trg->opt); //set pointer to reference target values
                reverse_complement(read_rc, read[LINE2], REVCOMP); //revcomp read sequence
                
                if ((anchr_pt = strstr(read_rc, SC1_anchr)) != NULL) {  //set SC1 anchor point
                    
                    end5p = anchr_pt + anchr_len; //the target RNA 5' end is the nt that follows the SC1 anchor seq
                    
                    //compare read sequence to target to assess whether the read is a true match
                    //the structure of this comparison considers any read that contains the target
                    //sequence as a substring that begins at index 0 to have mapped correctly, even
                    //if the read contains additional sequence downstream of the target substring.
                    //this allows users to ignore non-functional spacer sequence that is downstream
                    //of the target RNA during mapping by omitting it from the target sequences.
                    
                    for (rd_indx = 0, mtch = 1; end5p[rd_indx] && BC_val->tsq[rd_indx] && mtch; rd_indx++) {
                        if (end5p[rd_indx] != BC_val->tsq[rd_indx]) {
                            mtch = 0; //read is not a true match to the target
                        }
                    }
                    
                    //match fails if the read is shorter than the target,
                    //which is the case if the terminating null in end5p
                    //is reached before the terminating null in the target
                    if (!end5p[rd_indx] && BC_val->tsq[rd_indx]) {
                        mtch = 0;
                    }
                    
                    if (mtch) {
                        crnt_mpd_trg->cnt++;
                        
                        //determine if read is from bound or unbound channel
                        //the channel barcode is encoded in the first 4 bases of read 2. this information
                        //is removed from the read during fastp UMI processing and appended to the read id
                        //line (line 1) of the merged read
                        
                        channel = prcs_chnl_TDSPLY(&read[LINE1][0], met, &chnl_mtch_typ); //determine channel
                        met->chan_count[channel]++; //increment count for observed channel
                        
                        //record whether the mapped read originated in the bound or unbound channel
                        if (channel == BND) {
                            BC_val->chnl[BND]++; //increment bound counter for target
                        } else if (channel == UNB) {
                            BC_val->chnl[UNB]++; //increment unbound counter for target
                        } else {
                            BC_val->chnl[ERR]++; //increment error counter for target
                        }
                                                
                        if (testdata_MUX.run) { //if mapping test data, check that read mapped to correct target
                            compare_testdata_barcode_id(crnt_mpd_trg, &read[LINE1][0], &read[LINE2][0]);
                        }
                         
                    }
                }
            }
            //track ongoing processing. tracking is not done in testdata
            //mode because it interferes with error message printing
            if (i % 100000 == 0 && !testdata_MUX.run) {
                putchar('>');
                fflush(stdout);
            }
        }
    }
}

/* count_matched_compact_targets: count the number of copmact targets to which at least 1 read mapped */
int count_matched_compact_targets(compact_target * ctrg, target_params * trg_prms, int ctrg_cnt)
{
    int i = 0;               //general purpose index
    int matched_targets = 0; //number of targets with at least 1 mapped read
    
    for (i = 0; i < ctrg_cnt; i++) {         //for every target
        if (!ctrg[i].mul && ctrg[i].cnt) {   //if target is not redundant and has >=1 read mapped
            matched_targets++;               //increment matched targets count
        }
    }
    return matched_targets;  //return matched targets count
}

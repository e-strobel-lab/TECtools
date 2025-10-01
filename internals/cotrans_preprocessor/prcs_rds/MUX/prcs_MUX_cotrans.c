//
//  process_MUX_cotrans.c
//  
//
//  Created by Eric Strobel on 6/25/25.
//

#include <stdio.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../../variant_maker/make_barcodes.h"

#include "../../../TECdisplay_mapper/TECdisplay_mapper_defs.h"
#include "../../../TECdisplay_mapper/TECdisplay_mapper_structs.h"
#include "../../../TECdisplay_mapper/map_reads/map_expected/parse_mx_trgts.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"

#include "../UNV/call_fastp_TPROBE.h"
#include "../UNV/bypass_fastp.h"
#include "../UNV/prcs_chnl_TPROBE.h"
#include "../UNV/print_splitting_metrics.h"

#include "../../MUX_trgt_gen/mk_MUX_trgts.h"
#include "../../MUX_trgt_gen/mk_MUX_testdata.h"

#include "../MLT/prcs_MLT_cotrans.h"

#include "./testdataMUX_analysis.h"

#include "prcs_MUX_cotrans.h"

/* prcs_MUX_cotrans: manages processing of TECprobe-MUX data */
int prcs_MUX_cotrans(TPROBE_names * nm, FILE * fp_MUXtrgs, fastp_params fastp_prms, testdata_MUX_vars * testdata_MUX, int run_bypass_fastp)
{
    printf("processing multiplex cotrans\n");
        
    FILE *ifp[READ_MAX] = {NULL};  //pointers for input fastq files
    metrics met = {0};             //read processing metrics storage
    
    compact_target * ctrg  = NULL; //pointer to compact target structures
    opt_BC * BC_val = NULL;        //pointer to optional target value structures
    
    target_params trg_prms = {0};  //structure for storing target parameters
    TDSPLY_fasta wt = {0};         //storage for wt sequence information
    
    int brcd_cnt = 0;              //number of barcodes
    int ctrg_cnt = 0;              //number of compact targets stored
    int clcd_ctrg_cnt = 0;         //calculated number of barcode targets
    
    char line[MAX_LINE+1] = {0};   //array to store line
    
    /*parse_header_lines(fp_MUXtrgs, &trg_prms, &wt);
    printf("%d expected targets\n", trg_prms.xpctd);
    printf("%s: %s\n", wt.nm, wt.sq);
    abort();*/
    
    //count number of targets
    while (get_line(line, fp_MUXtrgs)) { //until all lines have been read
        if (line[0] == '>') {            //if reading first line of fasta entry
            get_line(line, fp_MUXtrgs);  //get the second line of the fasta entry
            brcd_cnt++;                  //increment barcode count
        } else {
            printf("prcs_MUX_cotrans: error - unexpected format for targets file. aborting...\n");
            abort();
        }
    }
    fclose(fp_MUXtrgs);                  //close targets file
    get_file(&(fp_MUXtrgs), nm->trgts);  //re-open targets file
    
    clcd_ctrg_cnt = brcd_cnt * 129;  //caculate expected number of barcode targets, including single subs and indels
    
    //allocate memory for compact targets
    if ((ctrg = calloc(clcd_ctrg_cnt, sizeof(*ctrg))) == NULL) {
        printf("prcs_MUX_cotrans: error - compact target memory allocation failed. aborting...\n");
        abort();
    }
    
    //allocate memory for barcode target optional values
    if ((BC_val = calloc(clcd_ctrg_cnt, sizeof(*BC_val))) == NULL) {
        printf("prcs_MUX_cotrans: error - optional target value memory allocation failed. aborting...\n");
        abort();
    }
    
    ctrg_cnt = mk_MUX_trgts(nm, ctrg, BC_val, fp_MUXtrgs, brcd_cnt, clcd_ctrg_cnt); //generate barcode targets
    met.srcTrgs = brcd_cnt; //record number of source barcodes
    met.targets = ctrg_cnt; //record number of targets generated
    
    /********* hash table initialization and construction **********/
    compact_h_node **htbl_MUX = NULL;                  //hash table root
    compact_h_node_bank hn_MUX_bank = {NULL, NULL, 0}; //bank for hash table nodes
    
    hn_MUX_bank.count = 0; //initialize hash table node count to zero
    
    //allocate TABLE_SIZE hash table node pointers
    if ((htbl_MUX = calloc(TABLE_SIZE, sizeof(*htbl_MUX))) == NULL) {
        printf("main: error - hash table memory allocation failed\n");
        abort();
    }
    
    //allocate BLOCK_SIZE hash table nodes
    if ((hn_MUX_bank.chn = calloc(BLOCK_SIZE, sizeof(*(hn_MUX_bank.chn)))) == NULL) {
        printf("main: error - hash table node memory allocation failed\n");
        abort();
    }
    
    compact_h_node_bank *crrnt_hn_MUX_bank = &hn_MUX_bank;          //pointer for handling hash node bank
    mk_htbl_MUX(htbl_MUX, crrnt_hn_MUX_bank, ctrg, ctrg_cnt, &met); //generate barcode target hash table
    /****** end of hash table initialization and construction ******/
    
    /*************** testdata generation *****************/
    if (testdata_MUX->run) {
        mk_MUX_testdata(nm, ctrg, ctrg_cnt); //generate test data
    }
    /*********** end of test data generation *************/
    
    /***************** obtain sample name prefix *******************/
    get_sample_name(nm->file[READ1], nm->smpl[READ1]); //get read 1 sample name
    get_sample_name(nm->file[READ2], nm->smpl[READ2]); //get read 2 sample name
    /************** end of obtain sample name prefix ***************/
    
    
    /************* process and split sequencing reads **************/
    mk_out_dir("split"); //make directory for output files
    
    if (testdata_MUX->run && run_bypass_fastp) {
        bypass_fastp(nm->file[READ1], nm->file[READ2], &ifp[0]);
    } else {
        call_fastp_TPROBE(nm->file[READ1], nm->file[READ2], &ifp[0], fastp_prms); //fastp pre-processing
    }
    
    //split input fastq by channel and barcode
    split_MUX_reads(&ifp[0], htbl_MUX, nm, ctrg, brcd_cnt, ctrg_cnt, &met, fastp_prms.mode);
    /********** end of process and split sequencing reads **********/
    
    print_splitting_metrics(nm, &met, fastp_prms); //print metrics
    if (testdata_MUX->run) {
        print_MUX_testdata_analysis(&met, &ctrg[0]);
    }
    
    system("rm ./split/R*out.fq"); //remove fastp output files
    system("gzip ./split/*.fq");   //compress split fastq files
    
    mk_config(nm, NULL, fastp_prms.mode); //make config file for run script generation
    
    return 1;
}

/* mk_htbl_MUX: makes compact target hash table */
/* hash table has linked list buckets for possible collisions */
void mk_htbl_MUX(compact_h_node ** htbl_MUX, compact_h_node_bank * bank, compact_target * ctrg, int count, metrics * met)
{
    extern uint64_t mutcode_mask; //mask to used to isolate barcode mutation code
    
    compact_h_node **p_rdnd = NULL; //pointer for h_node handling
    
    int i = 0;              //general purpose index
    int j = 0;              //general purpose index
    int new_node = 0;       //counts number of nodes that are assigned a target
    int redundant = 0;      //counts number of redundant targets (same seq as prev target)
    int collisions = 0;     //counts number of collisions between barcodes
    int ntv_collisions = 0; //counts number of native barcode collisions
    int blacklisted = 0;    //counts number of blacklisted mutant barcodes
    
    int error = 0; //flag that error occurred and program should be aborted
    
    uint64_t hash = 0;  //hash key
    
    char csq1[MAX_LINE+1] = {0};  //array for storing character-encoded barcode sequence
    char csq2[MAX_LINE+1] = {0};  //array for storing character-encoded barcode sequence
    
    opt_BC * p_BC_val_1 = NULL; //pointer to barcode target optional values
    opt_BC * p_BC_val_2 = NULL; //pointer to barcode target optional values
    
    char * type_str1 = NULL; //pointer to barcode type string for sequence 1
    char * type_str2 = NULL; //pointer to barcode type string for sequence 2
    int type_val1 = -1;      //sequence 1 barcode type value
    int type_val2 = -1;      //sequence 2 barcode type value
    
    for (i = 0; i < count; i++) {  //perform loop for each 3' end target
        
        hash = hash_brcd_trgt(&ctrg[i].bsq);                      //hash binary encoded sequence
        p_rdnd = srch_ctrg_htbl(&ctrg[i].bsq, hash, htbl_MUX, 0); //search hash table
        
        if ((*p_rdnd) == NULL) {                    //no existing hash table node for target sequence
            (*p_rdnd) = &bank->chn[bank->count++];  //assign node from hash node bank
            if (bank->count == BLOCK_SIZE) {        //check that bank was not filled
                extend_ch_bank(bank);               //extend bank if needed
            }
            (*p_rdnd)->ctrg = &(ctrg[i]);           //set node to point to current target
            new_node++;                             //increment new_node counter
            
        } else {
            ctrg[i].mul = 1; //set flag that current target is redudant with previous target
            redundant++;     //increment redundant counter
            
            //check whether redundant barcode mutants are linked to different reference barcodes
            if (check_brcd_diff((*p_rdnd)->ctrg, &ctrg[i])) {     //if brcds map to different ref brcds
                collisions++;                                     //increment collision counter
                
                bin2seq(csq1, &(*p_rdnd)->ctrg->bsq, MAX_LINE+1); //convert barcode 1 to string
                bin2seq(csq2, &ctrg[i].bsq, MAX_LINE+1);          //convert barcode 2 to string
                
                p_BC_val_1 = (opt_BC *)(*p_rdnd)->ctrg->opt; //set pointer to barcode 1 optional values
                p_BC_val_2 = (opt_BC *)ctrg[i].opt;          //set pointer to barcode 2 optional values
                
                type_str1 = get_target_type((*p_rdnd)->ctrg->bid & mutcode_mask, &type_val1);
                type_str2 = get_target_type(ctrg[i].bid & mutcode_mask, &type_val2);
                
                //report collision
                printf("\nCOLLISION:\n");
                
                if (!strcmp(csq1, csq2)) { //collided barcode sequences match (sanity check)
                    //printf("id1:  ");
                    //printbin((*p_rdnd)->ctrg->bid);
                    //printf("id2:  ");
                    //printbin(ctrg[i].bid);
                    printf("ref1: %s\nref2: %s\n", p_BC_val_1->ref->csq, p_BC_val_2->ref->csq);
                    printf("trg1: %s (%s)\ntrg2: %s (%s)\n\n", csq1, type_str1, csq2, type_str2);
                    
                } else { //collided barcode sequence do not match - this shouldn't be possible
                    printf("mk_htbl_mux: error - collision of non-identical barcodes. this should not be possible.\n");
                    printf("brcd: %s\n", csq1);
                    printf("brcd: %s\n", csq2);
                    printf("aborting...\n");
                    abort();
                }
                
                
                //TODO: since native barcodes are first and all native barcodes are guaranteed to be sufficiently distant from each other, I think only (*p_rdnd)->ctrg->bid needs to be tested here. Is there a case where ctrg[i].bid could be a native barcode?
                //test whether collision involved reference barcode
                if (!((*p_rdnd)->ctrg->bid & mutcode_mask) || !(ctrg[i].bid & mutcode_mask)) { //if either brcd was NAT

                    ntv_collisions++; //increment number of native barcode collisions
                    
                    //it should not be possible for a collision in which the second barcode was native to occur.
                    //native barcodes are added to the hash table first, and all native barcodes should be distinct
                    if (!(ctrg[i].bid & mutcode_mask)) { //if second barcode was native, throw error and abort
                        printf("\nmk_htbl_MUX: error - native collision in which second target was native barcode. this should not be possible. aborting after hash table generation...\n");
                        error = 1;
                    } else {                //otherwise
                        ctrg[i].bl++;       //blacklist second barcode
                        met->blacklisted++; //increment number of blacklisted barcodes
                        //note mul does not need to be set here since it is set above for all redundant
                        //targets regardless of whether a collision occurred
                    }
                    
                } else { //both barcode targets were mutants, set prev target as redundant, blacklist both
                    (*p_rdnd)->ctrg->mul = 1; //set previous target as redundant
                    (*p_rdnd)->ctrg->bl++;    //blacklist previous target
                    ctrg[i].bl++;             //blacklist current target
                    //as above, mul does not need to be set here since it is set above for all redundant
                    //targets regardless of whether a collision occurred
                }
            }
        }
    }
    
    //print hash table generation details
    printf("\n********************  hash table generation  ********************\n");
    printf("%6d barcode targets were assessed\n", i);
    printf("%6d barcode target sequences were assigned a node\n", new_node);
    printf("%6d redundant barcode target sequences were not assigned a node\n\n\n", redundant);
    
    if (error) { //if a serious error occurred above, throw error and abort
        printf("aborting due to fatal error. see message above for cause.\n");
        abort();
    }
    
    //prompt user regarding whether analysis should proceed
    //TODO: can likely remove this. if a collision occurs with a native barcode the program is stopped. collisions between mutant barcodes are unlikely to impact the results of the analysis since reads that map to mutant barcodes should be rare. instead, consider generating a file that reports all collisions so that there is a record
    char usr_resp[4] = {0};         //storage for user response
    char discard[MAX_LINE+1] = {0}; //array for flushing stdin
    int resp_provided = 0;          //flag that valid reponse was provided
    
    if (collisions) {
        printf("%d collision(s) between native/substitution/insertion/deletion barcode sequences were detected. proceed?\n", collisions);
        
        while (!resp_provided) {
            scanf("%3s", usr_resp);
            
            if (!strcmp(usr_resp, "yes")) {
                printf("proceeding with variant demultiplexing\n\n");
                resp_provided = 1;
                
            } else if (!strcmp(usr_resp, "no")) {
                printf("aborting variant demultiplexing\n\n");
                abort();
                
            } else {
                printf("invalid response. please enter \"yes\" or \"no\".\n");
                get_line(discard, stdin);
                
            }
        }
    }
}

/* hash_brcd_trgt: generates hash key for binary encoded sequence (up to 32 nt/64 bits */
uint64_t hash_brcd_trgt(binary_seq * bsq)
{
    return *bsq->sq % TABLE_SIZE;
}

/* check_brcd_diff: checks whether redundant barcode targets are linked to different barcode variants*/
int check_brcd_diff(compact_target * old, compact_target * new)
{
    if ((new->bid >> MUTCODE_BITS) != (old->bid >> MUTCODE_BITS)) {
        return 1;
    } else {
        return 0;
    }
}
/* map_brcd: map barcode to target using hash table */
compact_target * map_brcd(char * brcd_str, char * rc_brcd_str, compact_h_node **htbl_MUX, compact_target ** mpd_trg, metrics * met)
{
    extern struct testdata_MUX_vars testdata_MUX; //structure containing test data read analysis variables
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    compact_h_node ** p_rdnd = NULL; //pointer to hash table node pointer
    compact_target * ref_trg = NULL; //pointer to reference compact target structure
    
    opt_BC * BC_val  = NULL; //pointer to current barcode optional target value structure
    opt_BC * ref_val = NULL; //pointer to reference barcode optional target value structure
    
    binary_seq bsq = {0}; //variable for storing binary-encoded barcode sequence
    if ((bsq.sq = calloc(1, sizeof(*(bsq.sq)))) == NULL) {
        printf("map_brcd: error - failed to allocate memory for binary-encoded barcode sequence. aborting\n");
        abort();
    }
    
    uint64_t hash = 0; //hash value
    
    seq2bin_long(rc_brcd_str, &bsq, 1); //generate binary-encoded barcode sequence
    hash = hash_brcd_trgt(&bsq);        //hash binary-encoded barcode sequence
    
    //TODO: for now, trace search is off. add as debug option later
    p_rdnd = srch_ctrg_htbl(&bsq, hash, htbl_MUX, 0); //search hash table for match
    
    if ((*p_rdnd) != NULL) {                              //if found barcode match
        met->hits++;                                      //increment number of hits
        if (!strcmp(rc_brcd_str, (*p_rdnd)->ctrg->csq)) { //sanity check that barcode sequences match
            met->matches++;                               //if match, increment number of barcode sequence matches
        } else {                                          //otherwise, throw error and abort
            printf("map_brcd: error query string does not match has table entry. aborting...\n");
            abort();
        }
        
        if ((*p_rdnd)->ctrg->bl) { //if target is blacklisted,
            met->bl_mapped++;      //track mapped blacklisted read count
            return NULL;           //then return NULL
            
        } else { //otherwise track mapped read and return pointer to reference target
            met->mapped++;                           //track mapped read
            *mpd_trg = (*p_rdnd)->ctrg;              //set mapped target pointer to the current node's target
            
            BC_val = (opt_BC *)(*p_rdnd)->ctrg->opt; //set pointer to current target optional values
            BC_val->cnt++;                           //increment number of reads that mapped to current target
            switch (BC_val->typ) {                   //track number of NAT/SUB/INS/DEL reads
                case NAT: met->nat_cnt++; break;
                case SUB: met->sub_cnt++; break;
                case INS: met->ins_cnt++; break;
                case DEL: met->del_cnt++; break;
                default:
                    printf("map_brcd: unrecognized barcode mutant type. aborting...\n");
                    abort();
                    break;
            }
        
            ref_trg = BC_val->ref;            //set pointer to reference target of current mapped target
            ref_val = (opt_BC *)ref_trg->opt; //set pointer to reference target optional values
            ref_val->mpd++;                   //increment number of reads that mapped to the reference target
            free(bsq.sq);                     //free memory
            return ref_trg; //return reference target pointer
        }
    } else {             //otherwise
        met->unmapped++; //increment number of unmapped reads
        free(bsq.sq);    //free memory
        return NULL;     //return NULL
    }
}

/* split_MUX_reads: demultiplex TECprobe-MUX reads into separate fastq file */
void split_MUX_reads(FILE **ifp, compact_h_node **htbl_MUX, TPROBE_names * nm, compact_target * ctrg, int brcd_cnt, int ctrg_cnt, metrics * met, int mode)
{
    extern int debug;                             //flag to turn on debug mode
    extern struct testdata_MUX_vars testdata_MUX; //structure containing test data read analysis variables
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
        
    char out_nm[READ_MAX][(MAX_LINE*2)] = {{0}}; //arrays for output names
    
    int chnls2make = 2; //number of channels to make when splitting sequencing reads
    
    opt_BC * BC_val = NULL; //pointer to barcode target optional values
    
    //array for appending channel code to output fastq file names
    static const char chnl_code[CHANNEL_MAX][4] = {"UNT", "MOD", "ERR"};
    //UNT 0
    //MOD 1
    //ERR 2
    
    //generate output files
    for (i = 0; i < chnls2make; i++) {
        for (j = 0; j < brcd_cnt && j < MAX_BRCD_CNT; j++) {
            
            //format of output file names:
            //<sample name>_<barcode id>_<channel code>_R1/2.fq
            //<sample name> is obtained from the input file names and is read specific
            //<channel code> is UNT (untreated), MOD (modified), or ERR (error)
            //<barcode id> is the numerical id of the barcode
            sprintf(out_nm[READ1], "./split/%s_%s_%05llu_R1.fq", nm->smpl[READ1], chnl_code[i], (long long unsigned int)(ctrg[j].bid >> MUTCODE_BITS));
            sprintf(out_nm[READ2], "./split/%s_%s_%05llu_R2.fq", nm->smpl[READ2], chnl_code[i], (long long unsigned int)(ctrg[j].bid >> MUTCODE_BITS));
            
            BC_val = (opt_BC *)ctrg[j].opt; //set pointer to barcode optional values
            
            //open output files
            if ((BC_val->ofp[i][READ1] = fopen(out_nm[READ1], "w")) == NULL) {
                printf("split_MUX_readc: ERROR - could not generate %s file. try increasing ulimit -n to >(barcode count * 4). Aborting program...\n", out_nm[READ1]);
                abort();
            }
            if ((BC_val->ofp[i][READ2] = fopen(out_nm[READ2], "w")) == NULL) {
                printf("split_MUX_reads: ERROR - could not generate %s file. try increasing ulimit -n to >(barcode count * 4). Aborting program...\n", out_nm[READ2]);
                abort();
            }
        }
    }
    
    //general variables
    int got_line[READ_MAX] = {1};   //flag to indicate success of the get_line function
    int proceed = 1;                //flag to indicate that read processing should proceed
    
    //read variables
    char read1[FQ_LINES][MAX_LINE] = {{0}};  //array for storing all four lines of one read 1 fastq entry
    char read2[FQ_LINES][MAX_LINE] = {{0}};  //array for storing all four lines of one read 2 fastq entry
    char ipt_ID[READ_MAX][MAX_LINE] = {{0}}; //array for storing input lines when debug mode is on
    char ipt[READ_MAX][MAX_LINE] = {{0}};    //array for storing input reads for error checking
    
    //channel variables
    int  channel = -1;      //stores channel barcode code
    char chan_str[8] = {0}; //channel barcode string for appending attribute to read id
    
    //barcode variables
    char brcd_str[MAX_BARCODE_LEN+1] = {0};    //barcode string
    char rc_brcd_str[MAX_BARCODE_LEN+1] = {0}; //reverse complement of barcode string
    
    int r1_printed = -1; //flag to indicate that read 1 printing was successful
    int r2_printed = -1; //flag to indicate that read 2 printing was successful
    
    compact_target * crnt_ref_trg = NULL; //pointer to reference target of mapped target
    compact_target * crnt_mpd_trg = NULL; //pointer to mapped target
    
    for (i = 0; proceed; i++) {
        
        brcd_str[0] = chan_str[0] = '\0';   //initialize barcode and channel strings to 0
        channel = -1;                       //initialize channel value to -1
        crnt_ref_trg = crnt_mpd_trg = NULL; //initialize reference and mapped taret pointers to NULL
        
        for (j = 0; j < FQ_LINES; j++) {    //zero first index of fastq line storage
            read1[j][0] = read2[j][0] = '\0';
        }
        
        //copy fastq file lines for each read to read1 and read2 arrays
        for (j = 0; j < FQ_LINES && proceed; j++) {
            got_line[READ1] = get_line(&read1[j][0], ifp[READ1]); //get line for read 1
            got_line[READ2] = get_line(&read2[j][0], ifp[READ2]); //get line for read 2
            
            //copy read 1 and 2 sequences to ipt array for post-processing sequence verification
            if (j == LINE2 && got_line[READ1] && got_line[READ2]) {
                strcpy(ipt[READ1], read1[LINE2]); //store read 1 sequence for post-processing verification
                strcpy(ipt[READ2], read2[LINE2]); //store read 2 sequence for post-processing verification
            }
            
            //copy input ids for post-processing comparision when debug mode is on
            if (debug && j == LINE1 && got_line[READ1] && got_line[READ2]) {
                strcpy(ipt_ID[READ1], read1[LINE1]); //store read 1 ID for post-processing comparison
                strcpy(ipt_ID[READ2], read2[LINE1]); //store read 1 ID for post-processing comparison
            }
            
            if (!(got_line[READ1] && got_line[READ2])) { //test success of get_line for read 1 and read 2
                //test failed. this should only occur at end of file and get_line will
                //throw an error if an unexpected input line is encountered
                
                //test that the end of read 1 and 2 fastq files were reached at the same time
                if (got_line[READ1] == got_line[READ2]) { //reached EOF for both files
                    proceed = 0; //exit read processing loop
                } else {
                    printf("split_reads_3pEnd: error R1out.fq and R2out.fq (fastp output) do not contain the same number of lines. aborting...\n");
                    abort();
                }
            }
        }
        
        verify_read(read1); //verify read 1 integrity
        verify_read(read2); //verify read 2 integrity
        
        //determine channel and 3' end information for read
        if (proceed) {
            met->reads_processed++;                //track total number of reads processed
        
            //*** determine if read is from untreated or modified channel ***
            //the channel barcode is encoded in the first 5 bases of read 2. this information
            //is removed from the read during fastp UMI processing and appended to the read id
            //line (line 1) of reads 1 and 2
            
            channel = prcs_chnl_TPROBE(&read1[LINE1][0], met, mode);
            switch (channel) {
                case UNT: sprintf(chan_str, "UNT"); break;
                case MOD: sprintf(chan_str, "MOD"); break;
                case ERR: sprintf(chan_str, "ERR"); break;
                default:
                    printf("unexpected value for channel variable (%d). aborting...\n", channel);
                    abort();
                    break;
            }
            met->chan_count[channel]++; //increment count for observed channel
            
            get_brcd_str(brcd_str, &read1[LINE1][0]);                                     //get barcode from read1 id
            reverse_complement(rc_brcd_str, brcd_str, REVCOMP);                           //revcomp barcode string
            crnt_ref_trg = map_brcd(brcd_str, rc_brcd_str, htbl_MUX, &crnt_mpd_trg, met); //map barcode using hash table
            
            if (crnt_ref_trg != NULL) { //if barcode mapped
                
                if (testdata_MUX.run) { //if mapping test data, check that read mapped to correct target
                    compare_testdata_barcode_id(crnt_mpd_trg, &read2[LINE1][0], &read2[LINE2][0]);
                }
                
                BC_val = (opt_BC *)crnt_ref_trg->opt; //set pointer to ref barcode optional values
                
                // *** append attributes to read ids ***
                appnd_att(&read1[LINE1][0], &read2[LINE1][0], chan_str);          //append channel to read name
                appnd_att(&read1[LINE1][0], &read2[LINE1][0], crnt_ref_trg->cid); //append barcode id to read name
                if (debug) {printf("input:\t%s\n      \t%s\noutput:\t%s\n       \t%s\n\n------------\n\n",
                                   ipt_ID[READ1], ipt_ID[READ2],
                                   read1[LINE1], read2[LINE1]);}
                // *** end append attributes to read ids ***
                
                // *** verify read sequences ***

                //this is a sanity check, there is no reason
                //that the read sequences should ever change
                            
                if (!strcmp(ipt[READ1], read1[LINE2])) {
                    met->read_matches[READ1]++;
                } else {
                    printf("split_MUX_reads: CRITICAL - sequence integrity failure in read 1. this should NEVER happen! please contact estrobel@buffalo.edu regarding this error.\naborting...\n");
                    abort();
                }
                
                if (!strcmp(ipt[READ2], read2[LINE2])) {
                    met->read_matches[READ2]++;
                } else {
                    printf("split_MUX_reads: CRITICAL - sequence integrity failure in read 2. this should NEVER happen! please contact estrobel@buffalo.edu regarding this error.\naborting...\n");
                    abort();
                }
                // *** end verify read sequences ***
                
                if (channel > -1 && channel < ERR) { //if channel barcode corresponds to untreated or modified sample
                    
                    //print fastq entry to corresponding read 1/2 output files
                    for (j = 0; j < FQ_LINES; j++) {
                        r1_printed = r2_printed = -1; //set print success flags to -1
                        
                        r1_printed = fprintf(BC_val->ofp[channel][READ1], "%s\n", read1[j]);
                        r2_printed = fprintf(BC_val->ofp[channel][READ2], "%s\n", read2[j]);
                        
                        if (r1_printed <= 0 || r2_printed <= 0) {
                            printf("split_MUX_reads: error occurred when printing fastq file output. aborting...");
                            abort();
                        }
                    }
                } //TODO: consider adding option to keep discarded reads in separate file
                
                //track ongoing processing
                if (i+1 % 1000000 == 0) {
                    printf(">");
                }
            }
        }
    }

    printf("\n%d hits\n", met->hits);
    printf("%d matches\n\n", met->matches);
    
    printf("testdata bid match = %d\n", testdata_MUX.bid_match);
    
    printf("%d mapped\n", met->mapped);
    printf("%d unmapped\n\n", met->unmapped);
    
    printf("UNT TOT  = %d\n", met->chan_count[UNT]);
    printf("MOD TOT  = %d\n", met->chan_count[MOD]);
    printf("ERR TOT  = %d\n\n", met->chan_count[ERR]);
    printf("UNT FULL = %d\n", met->full_match[UNT]);
    printf("UNT PART = %d\n", met->part_match[UNT]);
    printf("MOD FULL = %d\n", met->full_match[MOD]);
    printf("MOD PART = %d\n", met->part_match[MOD]);
    
    
    //close output files
    for (i = 0; i < chnls2make; i++) {
        for (j = 0; j < brcd_cnt; j++) {
            
            BC_val = (opt_BC *)ctrg[j].opt;
            if ((fclose(BC_val->ofp[i][READ1])) == EOF || (fclose(BC_val->ofp[i][READ2])) == EOF) {
                printf("split_reads_3pEnd: error - error occurred when closing split fastq file. Aborting program...\n");
                abort();
            }
        }
    }
    
    return;
}

/* get_brcd_str: get barcode string from UMI in read ID */
void get_brcd_str(char * brcd_str, char * read1_ID)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    //iterate to first space in read 1 name, which indicates end of channel barcode
    for (i = 0; read1_ID[i] != ' ' && read1_ID[i]; i++) { ;}
    
    //check that read1_ID contained space and that expected UMI location is within array bounds
    if (!read1_ID[i]) { //check that loop did not exit on null character
        printf("prcs_chnl: error - unexpected read1 id line format. aborting...\n");
        abort();
    } else if (i <= ((MAX_BARCODE_LEN*2)+1)) { //check that negative index won't go outside of array bounds
        printf("prcs_chnl: error - unexpected short read1 id line. aborting...\n");
        abort();
    }
    
    //set index to start of barcode string, which is 2 barcode lengths + 1 (to account for an '_') upstream of the space
    i -= ((MAX_BARCODE_LEN*2)+1);
    
    if (read1_ID[i-1] != ':') { //if preceding character is not a colon, throw error and abort
        printf("get_brcd_str: error - unexpected format for trimmed read 1/2 head sequences. aborting...\n");
        abort();
    }
    
    //store barcode string
    for (j = 0; j < MAX_BARCODE_LEN && read1_ID[i] != '_'; i++, j++) {
        brcd_str[j] = read1_ID[i];
    }
    brcd_str[j] = '\0';
    
    if (read1_ID[i] != '_' || j != MAX_BARCODE_LEN) { //check that barcode was copied as expected
        printf("get_brcd_str: error - unexpected format for trimmed read 1/2 head sequences. aborting...\n");
        abort();
    }
    
    return;
}

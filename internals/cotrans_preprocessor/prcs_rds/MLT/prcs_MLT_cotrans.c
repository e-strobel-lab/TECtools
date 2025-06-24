//
//  prcs_MLT_cotrans.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/seq2bin_hash.h"
#include "parse_3pEnd_trgts.h"
#include "../UNV/call_fastp.h"
#include "../UNV/prcs_chnl.h"
#include "../UNV/print_splitting_metrics.h"
#include "../../../seq_utils/appnd_att.h"
#include "../../../seq_utils/isDNAbase.h"
#include "testdata3pEnd_analysis.h"
#include "mk_smooth_script.h"
#include "mk_MLT_config.h"
#include "../../../utils/debug.h"

#include "prcs_MLT_cotrans.h"

extern int debug;								  //flag to run debug mode
extern struct testdata_3pEnd_vars testdata_3pEnd; //structure containing test data read analysis variables

/*  prcs_MLT_cotrans: prcs_MLT_cotrans manage the functions that are required to split fastq files from multi-length cotranscriptional RNA structure probing experiments into separate fastq files based on the channel barcode and the RNA 3' end. */
int prcs_MLT_cotrans(TPROBE_names * nm, FILE * fp_3pEnd, fastp_params fastp_prms)
{
    FILE *ifp[READ_MAX] = {NULL};	//pointers for input fastq files
    metrics  met = {0};				//read processing metrics storage
    
    
    /***************** obtain sample name prefix *******************/
    get_sample_name(nm->file[READ1], nm->smpl[READ1]);
    get_sample_name(nm->file[READ2], nm->smpl[READ2]);
    /************** end of obtain sample name prefix ***************/
    
    
    /************* targets initialization and parsing **************/
    target *trgts = {NULL};			//pointer for array of hash table targets
    opt_3pEnd *end3p = {NULL};		//pointer for array of 3' end optional value structures
    target3p_params trg_prms = {0};	//structure for storing target parameters
    
    //end targets are only used for multi-length experiments
    //a targets array of BLOCK_SIZE (131072) should exceed the needs of normal use cases
    //TODO: consider adding an override to allow a larger targets array
    if (fastp_prms.mode == MULTI) {
        if ((trgts = calloc(BLOCK_SIZE, sizeof(*trgts))) == NULL) {
            printf("prcs_MLT_cotrans: error - hash table memory allocation failed\n");
            return 1;
        }
        if ((end3p = calloc(BLOCK_SIZE, sizeof(*end3p))) == NULL) {
            printf("prcs_MLT_cotrans: error - hash table memory allocation failed\n");
            return 1;
        }
        parse_3pEnd_trgts(fp_3pEnd, trgts, end3p, &trg_prms);
    } else if (fastp_prms.mode == SINGLE) {
        //set min and max target length to length of target
        trg_prms.min = nm->len;
        trg_prms.max = nm->len;
    }
    /********** end of targets initialization and parsing **********/

    
    /********* hash table initialization and construction **********/
    h_node **htbl_3end = NULL;			//hash table root
    h_node_bank hn_3end_bank = {NULL, NULL, 0};	//bank for hash table nodes
    
    hn_3end_bank.count = 0;
    if (fastp_prms.mode == MULTI) {
        //allocate TABLE_SIZE hash table node pointers
        if ((htbl_3end = calloc(TABLE_SIZE, sizeof(*htbl_3end))) == NULL) {
            printf("main: error - hash table memory allocation failed\n");
            return 1;
        }
        
        //allocate BLOCK_SIZE hash table nodes
        if ((hn_3end_bank.hn = calloc(BLOCK_SIZE, sizeof(*(hn_3end_bank.hn)))) == NULL) {
            printf("main: error - hash table node memory allocation failed\n");
            return 1;
        }
        
        h_node_bank *crrnt_hn_3end_bank = &hn_3end_bank; //pointer for handling hash node bank
        mk_htbl_3pEnd(htbl_3end, crrnt_hn_3end_bank, trgts, trg_prms.cnt); //generate 3' end target hash table
    }
    /****** end of hash table initialization and construction ******/
    
    
    /************* process and split sequencing reads **************/
    mk_out_dir("split");					//make directory for output files
    call_fastp(nm->file[READ1], nm->file[READ2], &ifp[0], fastp_prms);	//fastp pre-processing
    
    //split input fastq by channel and 3' end
    split_reads_3pEnd(&ifp[0], htbl_3end, nm, &met, trg_prms, fastp_prms.mode);
    /********** end of process and split sequencing reads **********/
    
    
    /********* print metrics and utility scripts, clean up *********/
    if (fastp_prms.mode == MULTI) {			//multi-length cotrans-specific processing
        mk_smooth_script(nm, trg_prms);		//generate smoothing script
        print_len_dist(&met, trg_prms);		//print the distribution of 3' end lengths
    }
    
    print_splitting_metrics(nm, &met, fastp_prms);	//print channel and 3' end pre-processing metrics
    if (testdata_3pEnd.run) {
        print_3pEnd_testdata_analysis(&met, trg_prms, trgts);	//print a report of test data read analysis
    }
    
    //close fastp output files
    if (fclose(ifp[READ1]) == EOF || fclose(ifp[READ2]) == EOF) {
        printf("split_reads_3pEnd: error - error occurred when closing fastp output files. Aborting program...\n");
        abort();
    }
    
    system("rm ./split/R*out.fq");			//remove fastp output files
    system("gzip ./split/*.fq");			//compress split fastq files
    
    mk_MLT_config(nm, trg_prms);			//generate config file for shapemapper 2 analysis
    
    /****** end of print metrics and utility scripts, clean up *****/
    
    return 1;
}



/* mk_htbl_3pEnd: construct hash table from target structure */
void mk_htbl_3pEnd(h_node **htbl, h_node_bank *bank, target *trgts, int count)
{
    /* hash table has linked list buckets for possible collisions */
    
    h_node **p_rdnd = NULL; //pointer for h_node handling
    
    int i = 0;
    int new_node = 0;		//counts number of nodes that are assigned a target
    int redundant = 0;		//counts number of redundant targets (same seq as prev target)
    
    for (i = 0; i < count; i++) {					//perform loop for each 3' end target
        p_rdnd = srch_htbl(trgts[i].key, htbl);		//search hash table for duplicate entries
        if ((*p_rdnd) == NULL) {					//no existing hash table node for target sequence
            (*p_rdnd) = &bank->hn[bank->count++];	//assign node from hash node bank
            if (bank->count == BLOCK_SIZE) {		//check that bank was not filled
                extend_h_bank(bank);				//extend bank if needed //TODO: currently not used
            }
            (*p_rdnd)->trg = &(trgts[i]);			//set node to point to current target
            new_node++;								//increment new_node counter
        } else {
            trgts[i].mul = 1;	//set flag that current target is redudant with previous target
            check_diff(*(*p_rdnd)->trg, trgts[i]);	//check length diff between query and existing target
            redundant++;							//increment redundant counter
        }
    }
    printf("\n********************  hash table generation  ********************\n");
    printf("%6d 3' end targets were assessed\n", i);
    printf("%6d 3' end target sequences were assigned a node\n", new_node);
    printf("%6d redundant 3' end target sequences were not assigned a node\n\n\n", redundant);
}



/* check_diff: check whether query target length matches existing target length.
 this is necessary because target mutations can result in the same
 target sequence being assigned to multiple transcript lengths.
 check_diff assesses whether the difference between the lengths of
 duplicate sequence entries are below the threshold set by TRG_DIF_THRESHOLD*/
void check_diff(target old, target new)
{
    int old_len = (*(opt_3pEnd *)old.opt).len; //dereference old length value
    int new_len = (*(opt_3pEnd *)new.opt).len; //dereference new length value

    if (abs(old_len - new_len) > TRG_DIF_THRESHOLD) {
        printf("mk_htbl_3pEnd: warning - non-unique targets with distance (%d) that exceeds TRG_DIF_THRESHOLD (%d) found in targets file\n", abs(old_len - new_len), TRG_DIF_THRESHOLD);
        printf("new:\t%s\t%s\nprev:\t%s\t%s\n", new.rc, new.id, old.rc, old.id);
    }
}



/* map_3pEnd: map RNA 3' end using hash table, return transcript length of the match.
 the sequence of the match is stored in end_trg_sq for use in clip_flanking */
int map_3pEnd(char * read, h_node **htbl, char * end_str, metrics * met, int trg_len)
{
    extern int debug;							//flag to turn on debug mode
    extern struct testdata_3pEnd_vars testdata_3pEnd;	//structure containing test data read analysis variables
    
    h_node **p_rdnd = NULL;						//pointer for hash table search result
    char rd_3pEnd[MAX_TRGT_SQ] = {0};			//stores read 1 head sequence to use for hash table search
    char rd_3pEnd_rc[MAX_TRGT_SQ] = {0};		//revcomp of R1 head, used for testdata warning
    
    //crnt_3pOpt is apointer to the opt_3pEnd struct that is currently being
    //assessed. this pointer simplifies the code be low since it reduces the
    //number of times that a void  pointer needs to be dereferenced.
    opt_3pEnd * crnt_3pOpt = NULL;
    
    //copy read1 head (RNA 3' end) to rd_3pEnd for hash table search
    int i = 0;
    for (i = 0; i < trg_len && read[i] && i < MAX_TRGT_SQ; i++) {
        rd_3pEnd[i] = read[i];
    }
    rd_3pEnd[i] = '\0';
    
    //if running debug mode, print report of 3' end mapping
    if (debug) {
        printf(">end determination\n");
        printf("read:\t%s\n", read);
        printf("query:\t%s\n", rd_3pEnd);
    }
    
    //search hash table for match to RNA 3' end
    p_rdnd = srch_htbl(rd_3pEnd, htbl);
    
    if ((*p_rdnd) != NULL) { //found RNA 3' end match
        
        (*p_rdnd)->trg->cnt++; //track the number of reads that map to each target
        
        //dereference optional void pointer in trg structure
        //to a opt_3pEnd structure pointer. this makes the
        //code that follows significantly more readable
        crnt_3pOpt = (opt_3pEnd *)(*p_rdnd)->trg->opt;
        
        //TODO: edit metrics output to report insertions and deletions too
        switch (crnt_3pOpt->typ) {
            case NAT: met->nat_cnt++; break;
            case SUB: met->sub_cnt++; break;
            case INS: met->ins_cnt++; break;
            case DEL: met->del_cnt++; break;
            default:
                printf("map_3pEnd: error - unrecognized 3' end mutant type. aborting...\n");
                abort();
                break;
        }
        
        //store transcript length as char string to append to read ids later
        sprintf(end_str, "%03d", crnt_3pOpt->len);
        
        //if running debug mode, print report of 3' end mapping
        if (debug) {
            printf("query:\t%s\t(repeated for comparison with node sequence)\n", rd_3pEnd);
            printf("node:\t%s\t", (*p_rdnd)->trg->rc);
            switch (crnt_3pOpt->typ) {
                case NAT: printf("NATIVE\n"); break;
                case SUB: printf("SUB\n"); break;
                case DEL: printf("DEL\n"); break;
                case INS: printf("INS\n"); break;
                default: break;
            }
            printf("node:\t%3d\n", crnt_3pOpt->len);
        }
        
        //sanity check that all mapped queries match the node sequence
        met->hits++;                                 //count number of hash table target hits
        if (!strcmp(rd_3pEnd, (*p_rdnd)->trg->rc)) { //count number of hits with expected seq (expect 100%)
            met->matches++;
        } else {
            printf("map_3pEnd: error - query string does not match hash table entry. aborting...\n");
            abort();
        }
        
        return crnt_3pOpt->len;
    } else { //no match for RNA 3' end in hash table
        
        //end string code of 000 indicates that the 3' end was unmappable
        sprintf(end_str, "000");
        if (debug) { printf("no match in hash table\n");}
        
        //if running test data analysis, print seed sequence information for unmapped reads
        if (testdata_3pEnd.run) {
            //
            //the reverse complement of rd_3pEnd is reported as the coding orientation of the
            //seed sequence ('sq=' below). this is because the head of read 1 is used to obtain
            //the seed sequence for 3' end mapping, but read 1 is the reverse complement of the
            //coding orientation of the target sequence.
            //
            //reporting this information is crucial because every test data read should be mappable
            //
            reverse_complement(rd_3pEnd_rc, rd_3pEnd, REVCOMP);
            printf("\nTEST DATA ANALYSIS: warning - no match in hash table for read\nseed: sq=%s\n      rc=%s\n",
                   rd_3pEnd_rc, rd_3pEnd);
        }
        return 0;
    }
}

/* verify_read: perform checks to assess the integrity of a read */
int verify_read(char (*rd)[MAX_LINE])
{
    int i = 0;
    
    int len_err = 0;
    int seq_err = 0;
    int qsc_err = 0;
    
    //verify that rd1 sequence length is equal to quality score length
    if (strlen(rd[LINE2]) != strlen(rd[LINE4])) {
        len_err = 1;
    }
    
    //verify that rd1 sequence and quality score lines contain only expected characters
    for (i = 0; rd[LINE2][i] && i < MAX_LINE; i++) {
        if (!isDNAbase(rd[LINE2][i]) && rd[LINE2][i] != 'N') {
            seq_err = 1;
        }
        
        if (rd[LINE4][i] < '!' || rd[LINE4][i] > 'J') {
            qsc_err = 1;
        }
    }
    
    if (len_err) {
        printf("verify_reads: error - read sequence and quality score lines are not the same length.\n");
        
    }
    if (seq_err) {
        printf("verify_reads: error - read sequence contains unexpected character %c.\n", rd[LINE2][i]);
        
    }
    if (qsc_err) {
        printf("verify_reads: error - read qscore contains unexpected character %c.\n", rd[LINE4][i]);
        
    }
    
    if (len_err || seq_err || qsc_err) {
        printf("\nread containing error(s):\n\n");
        for (i = 0; i < FQ_LINES; i++) {
            printf("%s\n", rd[i]);
        }
        printf("\naborting...\n");
        abort();
    }
    
    return 1;
}



/* split_reads_3pEnd: split reads from multi-length cotranscriptional RNA
 structure probing experiments into modified and untreated channels
 for each RNA 3' end length
 */
int split_reads_3pEnd(FILE **ifp, h_node **htbl, TPROBE_names * nm, metrics  * met, target3p_params trg_prms, int mode)
{
    printf("\nsplitting reads\n");
    
    extern int debug;							//flag to turn on debug mode
    extern struct testdata_3pEnd_vars testdata_3pEnd;	//structure containing test data read analysis variables
    
    int i = 0;
    int j = 0;
    
    FILE *out_fp[CHANNEL_MAX][END_MAX][READ_MAX] = {{{NULL}}};	//array for output file pointers
    char out_nm[READ_MAX][(MAX_LINE*2)] = {{0}}; 				//arrays for output names
    
    int chnls2make = 2;
    
    //array for appending channel code to output fastq file names
    static const char chnl_code[CHANNEL_MAX][4] = {"UNT", "MOD", "ERR"};
    //UNT 0
    //MOD 1
    //ERR 2
    
    //generate output files
    //
    //multi-length mode: generate UNT_R1, UNT_R2, MOD_R1, and MOD_R2 files for every
    //transcript length from min to max 3' end target. out_fp index corresponds directly
    //to transcript length for simplicity; this means that FILE pointers with an index
    //less than the minimum transcript length remain NULL.
    //
    //single-length mode: split input fastq into UNT_R1, UNT_R2, MOD_R1, and MOD_R2 files.
    //in single length mode these files are generated at the index for "length zero" and all
    //other FILE pointers remain NULL.
    for (i = 0; i < chnls2make; i++) {
        //in single length mode, this loop runs once because both trg_prms.min and
        //trg_prms.max are equal to zero
        for (j = trg_prms.min; j <= trg_prms.max; j++) {
            
            //format of output file names:
            //<sample name>_<channel code>_<transcript length>_R1/2.fq
            //<sample name> is obtained from the input file names and is read specific
            //<channel code> is UNT (untreated), MOD (modified), or ERR (error)
            //transcript length is a 3-digit value with leading zeros (multi) or '000' (single)
            sprintf(out_nm[READ1], "./split/%s_%s_%03d_R1.fq", nm->smpl[READ1], chnl_code[i], j);
            sprintf(out_nm[READ2], "./split/%s_%s_%03d_R2.fq", nm->smpl[READ2], chnl_code[i], j);
            
            if ((out_fp[i][j][READ1] = fopen(out_nm[READ1], "w")) == NULL) {
                printf("split_reads_3pEnd: ERROR - could not generate %s file. try increasing ulimit -n to >8000. Aborting program...\n", out_nm[READ1]);
                abort();
            }
            if ((out_fp[i][j][READ2] = fopen(out_nm[READ2], "w")) == NULL) {
                printf("split_reads_3pEnd: ERROR - could not generate %s file. try increasing ulimit -n to >8000. Aborting program...\n", out_nm[READ2]);
                abort();
            }
        }
    }
    
    //general variables
    int got_line[READ_MAX] = {1};	//flag to indicate success of the get_line function
    int proceed = 1;				//flag to indicate that read processing should proceed
    
    //read variables
    char read1[FQ_LINES][MAX_LINE] = {{0}};	 //array for storing all four lines of one read 1 fastq entry
    char read2[FQ_LINES][MAX_LINE] = {{0}};	 //array for storing all four lines of one read 2 fastq entry
    char ipt_ID[READ_MAX][MAX_LINE] = {{0}}; //array for storing input lines when debug mode is on
    char ipt[READ_MAX][MAX_LINE] = {{0}};	 //array for storing input reads for error checking
    
    //channel variables
    int  channel = -1;				//stores channel barcode code
    char chan_str[8] = {0};			//channel barcode string for appending attribute to read id
    
    //end variables
    int  end = 0;					//stores 3' end value
    char end_str[LEN_CODE] = {0};	//3' end value string for appending attribute to read id
    
    int r1_printed = -1; //flag to indicate that read 1 printing was successful
    int r2_printed = -1; //flag to indicate that read 2 printing was successful
    
    //process reads
    //
    //each iteration proceeds throught the following steps:
    //1. copy read 1 and read 2 fastq lines from the input files to the read1 and read2 arrays
    //2. determine the channel of the read
    //3. if multi-length, determine the 3' end of the read
    //4. append channel and transcript length info to read1 and read2 id lines (fastq line 1)
    //5. verify the integrity of read sequences
    //6. print the reads to an output fastq file based on channel and (if multi-length) 3' end
    //
    for (i = 0; proceed; i++) {
        
        end_str[0] = chan_str[0] = '\0'; //initialize end and channel strings to 0
        end = 0;		//initialize end value to 0
        channel = -1;	//initialize end value to -1
        
        for (j = 0; j < FQ_LINES; j++) {
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
            met->reads_processed++;				//track total number of reads processed
        
            /* determine if read is from untreated or modified channel */
            //the channel barcode is encoded in the first 5 bases of read 2. this information
            //is removed from the read during fastp UMI processing and appended to the read id
            //line (line 1) of reads 1 and 2
            
            channel = prcs_chnl(&read1[LINE1][0], met, mode);
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
            
            /* determine 3' end of transcript */
            if (mode == MULTI) {
                end = map_3pEnd(&read1[LINE2][0], htbl, end_str, met, trg_prms.sdLen);
                if (debug) { printf("return:\t%3d\nstring:\t%s\n\n", end, end_str);}
                if (end > 0) {
                    met->mapped++;		    //track number of mapped ends
                    met->len_dist[end]++;	//track distribution of mapped ends across transcripts
                } else {
                    met->unmapped++;	    //track number of unmapped ends
                }
                
                /* if running test data read analysis, track 3' end mapping accuracy */
                if (testdata_3pEnd.run) {
                    compare_testdata_3pEnd(end, read2[LINE1], read2[LINE2]);
                }
                
            } else if (mode == SINGLE) {
                end = nm->len; //set end to target length for use when printing output
                sprintf(end_str, "SGL"); //set end string code for single-length experiments
            }
            
            /****** append attributes to read ids ******/
            appnd_att(&read1[LINE1][0], &read2[LINE1][0], chan_str); //append channel to read name
            appnd_att(&read1[LINE1][0], &read2[LINE1][0], end_str);  //append 3' end  to read name
            if (debug) {printf("input:\t%s\n      \t%s\noutput:\t%s\n       \t%s\n\n------------\n\n",
                               ipt_ID[READ1], ipt_ID[READ2],
                               read1[LINE1], read2[LINE1]);}
            /**** end append attributes to read ids ****/
            
            /****** verify read sequences ******/

            //this is a sanity check, there is no reason
            //that the read sequences should ever change
                        
            if (!strcmp(ipt[READ1], read1[LINE2])) {
                met->read_matches[READ1]++;
            } else {
                printf("split_reads_3pEnd: CRITICAL - sequence integrity failure in read 1. this should NEVER happen! please contact estrobel@buffalo.edu regarding this error.\naborting...\n");
                abort();
            }
            
            if (!strcmp(ipt[READ2], read2[LINE2])) {
                met->read_matches[READ2]++;
            } else {
                printf("split_reads_3pEnd: CRITICAL - sequence integrity failure in read 2. this should NEVER happen! please contact estrobel@buffalo.edu regarding this error.\naborting...\n");
                abort();
            }
            /**** end verify read sequences ****/

            
            /* print reads to output fastq files
             mutli-length  mode requires a channel (UNT or MOD) and an end value >0
             single-length mode requires a channel (UNT or MOD) and that end=target length
             */

            if ((mode == MULTI  && channel > -1 && channel < ERR && end) ||
                (mode == SINGLE && channel > -1 && channel < ERR && end)) {
                
                for (j = 0; j < FQ_LINES; j++) {
                    r1_printed = r2_printed = -1; //set print success flags to -1
                    
                    r1_printed = fprintf(out_fp[channel][end][READ1], "%s\n", read1[j]);
                    r2_printed = fprintf(out_fp[channel][end][READ2], "%s\n", read2[j]);
                    
                    if (r1_printed <= 0 || r2_printed <= 0) {
                        printf("split_reads_3pEnd: error occurred when printing fastq file output. aborting...");
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
    printf("\n");
        
    //close output files
    for (i = 0; i < chnls2make; i++) {
        for (j = trg_prms.min; j <= trg_prms.max; j++) {
            if ((fclose(out_fp[i][j][READ1])) == EOF || (fclose(out_fp[i][j][READ2])) == EOF) {
                printf("split_reads_3pEnd: error - error occurred when closing split fastq file. Aborting program...\n");
                abort();
            }
        }
    }
    
    return 0;
}





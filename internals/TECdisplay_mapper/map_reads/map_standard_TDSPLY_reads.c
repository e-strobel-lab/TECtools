//
//  map_standard_TDSPLY_reads.c
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../TECdisplay_mapper_defs.h"
#include "../TECdisplay_mapper_structs.h"

#include "../../utils/gen_utils.h"
#include "../../utils/io_management.h"
#include "../../utils/debug.h"
#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"
#include "../../seq_utils/revcomp.h"
#include "../../seq_utils/basemap.h"
#include "../../seq_utils/mapping_metrics.h"
#include "./UNV/call_fastp_TDSPLY.h"
#include "./UNV/prcs_chnl_TDSPLY.h"
#include "./map_expected/get_key.h"
#include "./map_expected/parse_vmt_trgts.h"
#include "./map_expected/mk_htbl_TDSPLY.h"
#include "./map_expected/mk_chtbl_TDSPLY.h"
#include "./map_expected/mk_output_files.h"
#include "./map_expected/print_navigator_template.h"
#include "../testdata_analysis/mk_TDSPLY_test_data.h"
#include "../testdata_analysis/assess_TDSPLY_test_data.h"
#include "../testdata_analysis/parse_trgts_4_td_gen.h"

#include "./map_standard_TDSPLY_reads.h"

extern int debug;

/* prcs_standard_TDSPLY_reads: coordinates targets parsing, fastp processing, read mapping, and output file gen */
int prcs_standard_TDSPLY_reads(TDSPLY_names * nm, int trgt_ftype, char * minQ, fastp_params fastp_prms, testdata_vars * testdata, int mode, int use_ctrg)
{
    if (trgt_ftype != VMT_FILE) {
        printf("map_standard_TDSPLY_reads: error - expected vmt file input. aborting...\n");
        abort();
    }
    
    FILE *ifp = NULL;           //pointers for merged fastq file
    FILE * fp_trgs = NULL;      //pointer for targets file
    
    int i = 0;                  //general purpose index
    
    mapping_metrics met = {0};  //read processing metrics storage
    
    init_chnl_mtrcs_mem(&met, TDSPLY_CHANNEL_MAX); //initialize channel tracking memory
    
    /************* targets initialization and parsing **************/    
    target *refs = NULL;           //pointer for array of reference targets
    opt_ref *ref_val = NULL;       //pointer for array of optional reference target structures
    
    target *trgts = NULL;          //pointer for array of hash table standard targets
    compact_target *ctrg = NULL;   //pointer for array of hash table compact targets
    opt_mx_trg *trg_val = NULL;    //pointer for array of optional target structures
    
    target_params trg_prms = {0};  //structure for storing target parameters
    TDSPLY_fasta wt = {0};         //storage for wt sequence information
    
    //the variables below are used for testdata generation when mapping reads using compact targets
    //this circumvents the need to have a second function for testdata generation
    target * td_refs = NULL;         //pointer for array of reference targets
    opt_ref * td_ref_val = NULL;     //pointer for array of optional reference target structures
    target * td_trgts = NULL;        //pointer for array of hash table standard targets
    opt_mx_trg *td_trg_val = NULL;   //pointer for array of optional target structures
    target_params td_trg_prms = {0}; //structure for storing target parameters
    TDSPLY_fasta td_wt = {0};        //storage for wt sequence information
    
    int ret = 0; //variable for storing snprintf return value
    
    //variables for determining max number of variable bases across all reference targets
    char line[MAX_LINE+1] = {0}; //line storage
    char * p_vbs = NULL;         //pointer to "VBS:" in reference target line
    int ref_cnt = 0;             //number of reference targets
    int vbs_cnt = 0;             //number of variable bases in current reference target
    int max_vbs = 0;             //maximum number of variable bases observed across all reference targets
    
    h_node **htbl = NULL;                  //standard hash table root
    h_node_bank hn_bank = {NULL, NULL, 0}; //bank for standard hash table nodes
    hn_bank.count = 0;                     //initialize hash table bank count to 0
    
    compact_h_node **chtbl = NULL;                  //compact target hash table root
    compact_h_node_bank chn_bank = {NULL, NULL, 0}; //bank for compact target hash table nodes
    chn_bank.count = 0;                             //initialize compact hash table node count to zero
    compact_h_node_bank *crrnt_hn_bank = &chn_bank; //pointer for handling hash node bank
    
    //determine maximum number of variable bases
    get_file(&fp_trgs, nm->trgs); //open targets file
    
    while (get_line(line, fp_trgs)) {                      //for every line in the targets file
        if (!memcmp(line, "REF", strlen("REF"))) {         //if the line is a reference line
            ref_cnt++;                                     //increment reference line count
            
            if ((p_vbs = strstr(line, "VBS:")) != NULL) {  //set pointer to "VBS:" ahead of variable bases string
                
                //count the number of variable bases in the variable bases string of the current reference target
                for (i = strlen("VBS:"), vbs_cnt = 0; p_vbs[i] && !isspace(p_vbs[i]) && i < MAX_LINE; i++) {
                    if (isIUPACbase(p_vbs[i]) && (p_vbs[i+1] == '_' || isspace(p_vbs[i+1]))) {
                        vbs_cnt++;
                    }
                }
                
                if (vbs_cnt > max_vbs) { //current ref has more variable bases that previous refs
                    max_vbs = vbs_cnt;   //set new max_vbs
                }
            } else {
                printf("map_standard_TDSPLY_reads: error - failed to locate variable base sequence in vmt file reference target line. aborting...\n");
                abort();
            }
        }
    }
    fclose(fp_trgs); //close targets file
    
    //parse targets file header lines to get expected number of variants,
    //then allocate targets memory, then process targets file body and set targets
    
    get_file(&fp_trgs, nm->trgs);                //open targets file
    parse_header_lines(fp_trgs, &trg_prms, &wt); //parse targets file header lines
    printf("\nProcessing targets file that contains %d targets\n\n", trg_prms.xpctd);
    
    //allocate memory for reference targets and targets
    if ((refs = calloc(MAXREF, sizeof(*refs))) == NULL) {
        printf("map_standard_TDSPLY_reads: error - reference target memory allocation failed\n");
        return 1;
    }
    
    if ((ref_val = calloc(MAXREF, sizeof(*ref_val))) == NULL) {
        printf("map_standard_TDSPLY_reads: error - reference target value memory allocation failed\n");
        return 1;
    }
    
    //allocate standard target memory
    if ((trgts = calloc(trg_prms.xpctd, sizeof(*trgts))) == NULL) {
        printf("map_standard_TDSPLY_reads: error - target memory allocation failed\n");
        return 1;
    }
    
    if ((trg_val = calloc(trg_prms.xpctd, sizeof(*trg_val))) == NULL) {
        printf("map_standard_TDSPLY_reads: error - target value memory allocation failed\n");
        return 1;
    }
    
    if (max_vbs > SEQ2BIN_MAX_KEY || use_ctrg) {
        
        use_ctrg = 1; //set use_ctrg to 1 for use in code below
        
        printf("\nprocessing data using compact targets...\n\n");
        
        //allocate compact target memory
        if ((ctrg = calloc(trg_prms.xpctd, sizeof(*ctrg))) == NULL) {
            printf("map_standard_TDSPLY_reads: error - compact target memory allocation failed\n");
            return 1;
        }
        
        //store targets as compact targets and close targets file
        parse_vmt_trgts(fp_trgs, trgt_ftype, refs, ref_val,  ctrg, trg_val, &trg_prms, &wt, TDSPLY, CMPCT_TRGT_STRUCT);
        fclose(fp_trgs); //close targets file
        
        //allocate TABLE_SIZE hash table node pointers
        if ((chtbl = calloc(TABLE_SIZE, sizeof(*chtbl))) == NULL) {
            printf("map_barcoded_TDSPLY_reads: error - hash table memory allocation failed. aborting...\n");
            abort();
        }
        
        //allocate BLOCK_SIZE hash table nodes
        if ((chn_bank.chn = calloc(trg_prms.xpctd, sizeof(*(chn_bank.chn)))) == NULL) {
            printf("map_barcoded_TDSPLY_reads: error - hash table node memory allocation failed. aborting...\n");
            abort();
        }
        
        mk_chtbl_TDSPLY(chtbl, &chn_bank, ctrg, refs, &trg_prms);
        
    } else {
        
        printf("\nprocessing data using standard targets...\n\n");
        
        //store targets as standard targets and close targets file
        parse_vmt_trgts(fp_trgs, trgt_ftype, refs, ref_val, trgts, trg_val, &trg_prms, &wt, TDSPLY, STD_TRGT_STRUCT);
        fclose(fp_trgs); //close targets file
        
        //allocate TABLE_SIZE hash table node pointers
        if ((htbl = calloc(TABLE_SIZE, sizeof(*htbl))) == NULL) {
            printf("main: error - hash table memory allocation failed\n");
            return 1;
        }
        
        //allocate trg_prms.xpctd hash table nodes
        if ((hn_bank.hn = calloc(trg_prms.xpctd, sizeof(*(hn_bank.hn)))) == NULL) {
            printf("main: error - hash table node memory allocation failed\n");
            return 1;
        }
        
        mk_htbl_TDSPLY(htbl, &hn_bank, trgts, refs, &trg_prms); //generate target hash table
    }
    
    /****** end of hash table initialization and construction ******/
    
    
    /*************** testdata generation *****************/
    if (mode == MAP_TEST_DATA) {
        if (use_ctrg) {
            parse_trgts_4_td_gen(nm, trgt_ftype, &td_refs, &td_ref_val, &td_trgts, &td_trg_val, &td_trg_prms, &td_wt);
            mk_TDSPLY_test_data(nm, td_refs, td_trgts, &trg_prms);
        } else {
            mk_TDSPLY_test_data(nm, refs, trgts, &trg_prms); //generate test data
        }
    }
    /*********** end of test data generation *************/
    
    /***************** obtain sample name prefix *******************/
    get_sample_name(nm->file[READ1], nm->smpl[READ1]); //get read 1 sample name
    get_sample_name(nm->file[READ2], nm->smpl[READ2]); //get read 2 sample name
    ret = snprintf(nm->mrg, MAX_LINE, "%s_merged", nm->smpl[READ1]);    //construct merged read sample name
    if (ret >= MAX_LINE || ret < 0) {
        printf("map_standard_TDSPLY_reads: error - error when constructing output file name. aborting...\n");
        abort();
    }
    /************** end of obtain sample name prefix ***************/
    
    /************* process sequencing reads **************/
    mk_out_dir("processed");           //make directory for output files
    call_fastp_TDSPLY(nm, fastp_prms); //fastp pre-processing todo pass by reference?
    
    char processed_file[MAX_LINE]; //array for storing fastp-processed file name
    char gunzip[MAX_LINE];         //array for storing gunzip command
    
    //construct merged read file name
    ret = snprintf(processed_file, MAX_LINE, "./processed/%s.fq", nm->mrg);
    if (ret >= MAX_LINE || ret < 0) {
        printf("map_standard_TDSPLY_reads: error - error when constructing merged read file name. aborting...\n");
        abort();
    }
    
    //construct gunzip command
    ret = snprintf(gunzip, MAX_LINE, "gunzip %s", processed_file);
    if (ret >= MAX_LINE || ret < 0) {
        printf("map_standard_TDSPLY_reads: error - error when constructing gunzip command. aborting...\n");
        abort();
    }
    
    system(gunzip); //gunzip merged reads file
    
    //open merged read file
    if ((ifp = fopen(processed_file, "r")) == NULL) { //open merged reads file
        printf("error: could not open %s as merged read file. Aborting program...\n", processed_file);
        abort();
    }
    
    //map reads
    if (use_ctrg) {
        map_standard_TDSPLY_reads(ifp, chtbl, refs, ctrg, minQ, &trg_prms, &met, testdata, mode, use_ctrg); //map reads to targets
        set_utility_target(ctrg, trgts, trg_prms.xpctd);
    } else {
        map_standard_TDSPLY_reads(ifp, htbl, refs, trgts, minQ, &trg_prms, &met, testdata, mode, use_ctrg); //map reads to targets
    }
    
    trg_prms.mapped2 = count_matched_targets(trgts, &trg_prms); //count non-redundant targets with >=1 mapped read
    
    print_output(trgts, &trg_prms, nm, fastp_prms.mode); //print output file
    print_navigator_template(refs, &wt, &trg_prms);      //print navigator template file
    
    //close merged read file
    if ((fclose(ifp)) == EOF) {
        printf("map_standard_TDSPLY_reads: error - error occurred when closing merged read file. Aborting program...\n");
        abort();
    }
    
    /********** end of process sequencing reads **********/
    
    /************ print metrics and clean up *************/
    print_metrics(&trg_prms, &met, nm); //print mapping metrics
    
    if (testdata->run) {                       //print test data analysis
        print_testdata_analysis(&met, testdata, &trg_prms, trgts);
    }
    
    //compress merged fastq file
    char gzip[MAX_LINE]; //array to store gzip command
    
    //construct gzip command
    ret = snprintf(gzip, MAX_LINE, "gzip %s", processed_file);
    if (ret >= MAX_LINE || ret < 0) {
        printf("map_standard_TDSPLY_reads: error - error when constructing gzip command. aborting...\n");
        abort();
    }
    system(gzip); //gzip merged read file
    /*********** end of print metrics clean up ***********/
    
    //free memory
    free(refs);
    free(trgts);
    free(ref_val);
    free(trg_val);
    free(htbl);
    free(hn_bank.hn);
    
    return 1;
}

/* map_standard_TDSPLY_reads: map reads to user-supplied vmt targets */
void map_standard_TDSPLY_reads(FILE *ifp, void *htbl, target *refs, void *trgts, char * minQ, target_params * trg_prms, mapping_metrics * met, testdata_vars * testdata, int mode, int use_ctrg)
{
    extern int debug;  //flag to turn on debug mode
    
    //anchor variables
    const char SC1_anchr[15] = "GGCCTTCGGGCCAA";   //structure cassette 1 sequence; anchor for finding RNA 5' end
    int anchr_len = strlen(SC1_anchr);             //length of SC1 sequence
        
    //general variables
    int i = 0;                             //general purpose index
    int j = 0;                             //general purpose index
    int proceed = 1;                       //flag to indicate that read processing should proceed
    
    //read variables
    char read[FQ_LINES][MAX_LINE] = {{0}}; //array for storing all four lines of fastq entry
    char read_rc[MAX_LINE] = {0};          //array for storing reverse complement of read
    char qscore_rev[MAX_LINE] = {0};       //array for storing reversed qscore string
    
	//mapping variables
    char * anchr_pt = NULL;                //pointer to the location of the SC1 hairpin anchor sequence
    char * end5p = NULL;                   //pointer to 5' end of the RNA target sequence (after SC1 hairpin)
    
    uint64_t offset = 0;                   //offset to start of end5p in read_rc, used to set qscore pointer
    char * qscore5p = NULL;                //pointer to qscores starting at 5' end of native sequence
        
    /* testdata variables */
    char id_line[MAX_LINE] = {0}; //array to store testdata id line (fastq line 1)
    char * td_trg_id = NULL;      //pointer to expected target id in test data id line
    int td_ref_indx = -1;         //test data reference index
    int crnt_mut_cd = -1;         //current mutation code
    
    /* process reads
     read are processed as follows:
     1. store current fastq entry in the read array
     2. revcomp read and reverse qscore so that the sequence is in the sense orientation
     3. locate target 5' end within read by searching for SC1 anchor sequence
     4. generate key based on variable base positions in the reference sequence and attempt to map read
     5. if read does not map and there is another reference sequence, try mapping with the next reference
     6. when a read maps, determine the channel from which the read originates
     7. increment the appropriate channel counter for the target to which the read mapped
     */
    for (i = 0; proceed; i++) {
        
        /*** zero critical values before each loop iteration ***/
        for (j = 0; j < FQ_LINES; j++) { //initialize fastq arrays to 0
            read[j][0] = '\0';
        }
        
        read_rc[0] = '\0';    //initialize read revcomp array to 0
        qscore_rev[0] = '\0'; //initialize qscore reverse array to 0
        
        anchr_pt = NULL;      //initialize anchor point pointer to NULL
        end5p = NULL;         //initialize 5' end pointer to NULL
        
        offset = 0;           //initialize qscore start offset to 0
        qscore5p = NULL;      //initialize qscore5p pointer to NULL
                
        id_line[0] = '\0';    //initialize testdata id line to 0
        td_trg_id = NULL;     //initialize target id pointer to NULL
        td_ref_indx = -1;     //initialize testdata reference index to -1
        crnt_mut_cd = -1;     //initialize current mutation code to -1
        
        /* copy fastq file lines to read array */
        for (j = 0; j < FQ_LINES && proceed; j++) {
            if (!get_line(&read[j][0], ifp)) { //get merged read fastq line, test success of get_line
                //test failed. this should only occur at end of file and get_line will
                //throw an error if an unexpected input line is encountered
                proceed = 0;
            }
        }
        
        /* map read */
        if (proceed) {
            met->reads_processed++; //track total number of reads processed
            
            if (testdata->run) {              //if mapping test data...
                strcpy(id_line, read[LINE1]); //make copy of id line
                td_ref_indx = parse_testdata_id(testdata, &td_trg_id, &crnt_mut_cd, id_line); //parse test data id line
            }
        
            /* reverse complement read and reverse qscore */
            //mapping is done with the reverse complement of the read. although revcomping slows
            //down the analysis, it is much simpler to identify the anchor point and generate
            //a hash key using the reverse complement of the read, which is the sense orientation
            //of the sequence.
            reverse_complement(read_rc, read[LINE2], REVCOMP);
            reverse_complement(qscore_rev, read[LINE4], REVERSE);
            
            /* map the anchor point and set the 5' end of the target RNA */
            //the SC1 anchor immediately precedes the start of the RNA of interest.
            if ((anchr_pt = strstr(read_rc, SC1_anchr)) != NULL) {	//set SC1 anchor point
                end5p = anchr_pt + anchr_len;   //the target RNA 5' end is the nt that follows the SC1 anchor sequence
                offset = (uint64_t)(&end5p[0]) - (uint64_t)(&read_rc[0]); //determine end5p offset from read_rc[0]
                qscore5p = &qscore_rev[offset]; //set pointer to start of target RNA qscore
                
                if (debug) {
                    printf(">target RNA substring identification\n");
                    //printf("merged read input:\n%s\n%s\n\n", read[LINE2], read[LINE4]);
                    //printf("merged read revcomp:\n%s\n%s\n\n", read_rc, qscore_rev);
                    printf("anchor:\t%s\n", anchr_pt);
                    printf("5p end:\t");
                    for (j = 0; j < anchr_len; j++) { putchar(' ');}
                    printf("%s\n", end5p);
                    printf("qscore:\t");
                    for (j = 0; j < anchr_len; j++) { putchar(' ');}
                    printf("%s\n\n", qscore5p);
                }
                
                if (use_ctrg) {
                    map_read_to_cmpct_trgt(end5p, qscore5p, (compact_h_node **)htbl, refs, (compact_target *)trgts, minQ, trg_prms, met, testdata, &read[LINE1][0], td_trg_id, td_ref_indx, crnt_mut_cd, mode);
                } else {
                    map_read_to_std_trgt(end5p, qscore5p, (h_node **)htbl, refs, (target *)trgts, minQ, trg_prms, met, testdata, &read[LINE1][0], td_trg_id, td_ref_indx, crnt_mut_cd, mode);
                }
            }
            
            //track ongoing processing. tracking is not done in testdata
            //mode because it interferes with error message printing
            if (i % 100000 == 0 && !testdata->run) {
                putchar('>');
                fflush(stdout);
            }
        }
    }
}

/* map_read_to_std_trgt: map reads to target using a key with a max char length of 32 */
void map_read_to_std_trgt(char * end5p, char * qscore5p, h_node **htbl, target *refs, target *trgts, char * minQ, target_params * trg_prms, mapping_metrics * met, testdata_vars * testdata, char * rd_id, char * td_trg_id, int td_ref_indx, int crnt_mut_cd,  int mode) {
    
    extern int debug;  //flag to turn on debug mode
    
    if (debug) {printf(">read mapping\n");}
    
    //map the read using the target hash table
    //targets can be generated from multiple reference sequences with different key structures.
    //therefore, reads are mapped by sequentially generating keys using each reference sequence
    //as a template until all reference sequences are tested or a match is found.
    
    int ref_indx = 0;                      //index for reference targets
    int rd_indx = 0;                       //index for end5p array
    int mtch = 0;                          //flag that read is a match to target
    
    char key[SEQ2BIN_MAX_KEY+1] = {0};     //array for storing variable base string for use as hash table key
    h_node **p_rdnd = NULL;                //pointer to hash table node
    opt_mx_trg * crnt_trg_val = NULL;      //pointer to TECdisplay-specific target values

    int  channel = -1;                     //stores channel barcode code
    int  chnl_mtch_typ = -1;               //stores channel match type (full/partial)
    
    for (ref_indx = 0, mtch = 0; ref_indx < trg_prms->r_cnt && ref_indx < MAXREF && !mtch; ref_indx++) {
                
        if (get_key(key, end5p, qscore5p, &minQ[Q_VARIABLE], &refs[ref_indx], READ_KEY, MAX32CHAR_TARGET)) { //generate key

            if (debug) {printf("ref%d key: %s\n", ref_indx, key);}
            
            p_rdnd = srch_htbl(key, htbl);  //search hash table for match to hash key
            if ((*p_rdnd) != NULL) {        //if hash key match is found, perform full sequence comparison
                met->hits++;                //track number of hash key hits
                if (debug) {printf("hash table hit\n");}
                
                //compare read sequence to target to assess whether the read is a true match
                //the structure of this comparison considers any read that contains the target
                //sequence as a substring that begins at index 0 to have mapped correctly, even
                //if the read contains additional sequence downstream of the target substring.
                //this allows users to ignore non-functional spacer sequence that is downstream
                //of the target RNA during mapping by omitting it from the target sequences.
                for (rd_indx = 0, mtch = 1; end5p[rd_indx] && (*p_rdnd)->trg->sq[rd_indx] && mtch; rd_indx++) {
                    if (end5p[rd_indx] != (*p_rdnd)->trg->sq[rd_indx]) {
                        mtch = 0; //read is not a true match to the target
                    }
                }
                
                //match fails if the read is shorter than the target,
                //which is the case if the terminating null in end5p
                //is reached before the terminating null in the target
                if (!end5p[rd_indx] && (*p_rdnd)->trg->sq[rd_indx]) {
                    mtch = 0;
                }
                
                //if the min constant base qscore is higher than 0 (!), check that
                //all bases in the target RNA segment of the read meet or exceed
                //the minimum constant base qscore
                if (minQ[Q_CONSTANT] > '!') {
                    if (!test_cbase_qscores(qscore5p, &minQ[Q_CONSTANT], &refs[ref_indx])) {
                        mtch = 0;
                    }
                }
                
                if (mtch) { //read is a true match to the target
        
                    (*p_rdnd)->trg->cnt++;                            //increment target match counter
                    crnt_trg_val = (opt_mx_trg *)(*p_rdnd)->trg->opt; //set pointer to target values
                    met->matches++;                                   //track number of true matches
                    if (debug) {printf("read mapped using ref%d key\n\n", ref_indx);}
                    
                    /* determine if read is from bound or unbound channel */
                    //the channel barcode is encoded in the first 4 bases of read 2. this information
                    //is removed from the read during fastp UMI processing and appended to the read id
                    //line (line 1) of the merged read
                    channel = prcs_chnl_TDSPLY(rd_id, met, &chnl_mtch_typ); //determine channel
                    met->chan_count[channel]++;                             //increment count for observed channel
                    
                    //record whether the mapped read originated in the bound or unbound channel
                    if (channel == BND) {
                        crnt_trg_val->bnd++; //increment bound counter for target
                    } else if (channel == UNB) {
                        crnt_trg_val->unb++; //increment unbound counter for target
                    }
                    
                    if (mode == MAP_TEST_DATA) { //if mapping test data
                        
                        //evaluate the testdata match. if a mutant read maps to a target that was
                        //generated from a different variant template (return = -1), run
                        //crrct_testdata_nonsrc_mtch. this decrements all match counter variables
                        //so that the mutant read is ignored when assessing expected mapping outcomes.
                        //this correction is necessary for running testdata analysis using targets files
                        //that were generated from multiple closely related sequences in which a mutant
                        //of one variant could map to a target of another variant
                        
                        if (eval_testdata_mtch(testdata, td_ref_indx, td_trg_id, crnt_mut_cd, end5p, refs, (opt_mx_trg *)((*p_rdnd)->trg->opt), (*p_rdnd)->trg->id, (*p_rdnd)->trg->sq)) {
                            crrct_aberrant_testdata_mtch(&((*p_rdnd)->trg->cnt), crnt_trg_val, met,channel, chnl_mtch_typ, testdata);
                        }
                    }
                } else {
                    if (debug) {printf("hash table hit was not a true match\n\n");}
                }
            } else {
                if (debug) {printf("no hash table hit\n\n");}
            }
        } else {
            if (debug) {printf("ref%d key: key generation failed\n\n", ref_indx);}
        }
    }
}

/* map_read_to_cmpct_trgt: map reads to compact target, which may contain a key that is longer than 32 chars */
void map_read_to_cmpct_trgt(char * end5p, char * qscore5p, compact_h_node **chtbl, target *refs, compact_target *ctrg, char * minQ, target_params * trg_prms, mapping_metrics * met, testdata_vars * testdata, char * rd_id, char * td_trg_id, int td_ref_indx, int crnt_mut_cd,  int mode) {
    
    extern int debug;  //flag to turn on debug mode
    
    if (debug) {printf(">read mapping\n");}
    
    //map the read using the target hash table
    //targets can be generated from multiple reference sequences with different key structures.
    //therefore, reads are mapped by sequentially generating keys using each reference sequence
    //as a template until all reference sequences are tested or a match is found.
    
    int ref_indx = 0;                      //index for reference targets
    int rd_indx = 0;                       //index for end5p array
    int mtch = 0;                          //flag that read is a match to target
    
    char key[MAX_LINE+1] = {0};            //array for storing variable base string for use as hash table key
    binary_seq query_bsq = {0};            //structure for storing binary-encoded query sequence
    uint64_t hash = 0;                     //hash table key generated from query_bsq
    
    compact_h_node **p_rdnd = NULL;        //pointer to hash table node
    opt_mx_trg * crnt_trg_val = NULL;      //pointer to TECdisplay-specific target values

    int  channel = -1;                     //stores channel barcode code
    int  chnl_mtch_typ = -1;               //stores channel match type (full/partial)
    
    for (ref_indx = 0, mtch = 0; ref_indx < trg_prms->r_cnt && ref_indx < MAXREF && !mtch; ref_indx++) {
                
        if (get_key(key, end5p, qscore5p, &minQ[Q_VARIABLE], &refs[ref_indx], READ_KEY, UNBOUNDED_TARGET)) { //generate key
            
            if (debug) {printf("ref%d key: %s\n", ref_indx, key);}
            
            seq2bin_long(key, &query_bsq, MAX_BLOCKS); //encode key as bin_seq
            hash = hash_long_bsq(&query_bsq);          //generate hash table key
            
            p_rdnd = srch_ctrg_htbl(&query_bsq, hash, chtbl, 0);  //search hash table for match to hash key
            
            if ((*p_rdnd) != NULL) {        //if hash key match is found, perform full sequence comparison
                met->hits++;                //track number of hash key hits
                if (debug) {printf("hash table hit\n");}
                
                //compare read sequence to target to assess whether the read is a true match
                //the structure of this comparison considers any read that contains the target
                //sequence as a substring that begins at index 0 to have mapped correctly, even
                //if the read contains additional sequence downstream of the target substring.
                //this allows users to ignore non-functional spacer sequence that is downstream
                //of the target RNA during mapping by omitting it from the target sequences.
                for (rd_indx = 0, mtch = 1; end5p[rd_indx] && (*p_rdnd)->ctrg->csq[rd_indx] && mtch; rd_indx++) {
                    if (end5p[rd_indx] != (*p_rdnd)->ctrg->csq[rd_indx]) {
                        mtch = 0; //read is not a true match to the target
                    }
                }
                
                //match fails if the read is shorter than the target,
                //which is the case if the terminating null in end5p
                //is reached before the terminating null in the target
                if (!end5p[rd_indx] && (*p_rdnd)->ctrg->csq[rd_indx]) {
                    mtch = 0;
                }
                
                //if the min constant base qscore is higher than 0 (!), check that
                //all bases in the target RNA segment of the read meet or exceed
                //the minimum constant base qscore
                if (minQ[Q_CONSTANT] > '!') {
                    if (!test_cbase_qscores(qscore5p, &minQ[Q_CONSTANT], &refs[ref_indx])) {
                        mtch = 0;
                    }
                }
                
                if (mtch) { //read is a true match to the target
        
                    (*p_rdnd)->ctrg->cnt++;                            //increment target match counter
                    crnt_trg_val = (opt_mx_trg *)(*p_rdnd)->ctrg->opt; //set pointer to target values
                    met->matches++;                                    //track number of true matches
                    if (debug) {printf("read mapped using ref%d key\n\n", ref_indx);}
                    
                    /* determine if read is from bound or unbound channel */
                    //the channel barcode is encoded in the first 4 bases of read 2. this information
                    //is removed from the read during fastp UMI processing and appended to the read id
                    //line (line 1) of the merged read
                    channel = prcs_chnl_TDSPLY(rd_id, met, &chnl_mtch_typ); //determine channel
                    met->chan_count[channel]++;                             //increment count for observed channel
                    
                    //record whether the mapped read originated in the bound or unbound channel
                    if (channel == BND) {
                        crnt_trg_val->bnd++; //increment bound counter for target
                    } else if (channel == UNB) {
                        crnt_trg_val->unb++; //increment unbound counter for target
                    }
                    
                    if (mode == MAP_TEST_DATA) { //if mapping test data
                        
                        //evaluate the testdata match. if a mutant read maps to a target that was
                        //generated from a different variant template (return = -1), run
                        //crrct_testdata_nonsrc_mtch. this decrements all match counter variables
                        //so that the mutant read is ignored when assessing expected mapping outcomes.
                        //this correction is necessary for running testdata analysis using targets files
                        //that were generated from multiple closely related sequences in which a mutant
                        //of one variant could map to a target of another variant
                        
                        if (eval_testdata_mtch(testdata, td_ref_indx, td_trg_id, crnt_mut_cd, end5p, refs, (opt_mx_trg *)((*p_rdnd)->ctrg->opt), (*p_rdnd)->ctrg->cid, (*p_rdnd)->ctrg->csq)) {
                            crrct_aberrant_testdata_mtch(&((*p_rdnd)->ctrg->cnt), crnt_trg_val, met, channel, chnl_mtch_typ, testdata);
                        }
                    }
                } else {
                    if (debug) {printf("hash table hit was not a true match\n\n");}
                }
            } else {
                if (debug) {printf("no hash table hit\n\n");}
            }
        } else {
            if (debug) {printf("ref%d key: key generation failed\n\n", ref_indx);}
        }
    }
}


/* test_cbase_qscores: test whether all constant bases meet or exceed the minimum qscore */
int test_cbase_qscores(char * qscore5p, char * minQc, target *refs)
{
    int i = 0;                  //general purpose index
    int j = 0;                  //general purpose index
    int len = strlen(refs->sq); //length of reference target
    int pass = 1;               //flag that all qscores passed, initalized to true
    
    opt_ref * crnt_ref_val = (opt_ref *)refs->opt; //pointer to reference target optional values
    
    //test whether constant bases meet/exceed the minimum qscore
    for (i = 0; i < len && qscore5p[i] && pass; i++) {
        
        //check whether the current index matches a vbase index at each position. the
        //vbase indices are in ascending order, so it is only necessary to check whether
        //the next vbase index has been reached yet
        
        if (i != crnt_ref_val->vb_pos[j]) {    //if the current base is a constant base
            
            if (qscore5p[i] < *minQc) {        //if min cbase qscore is not met/exceeded
                pass = 0;                      //set pass to false
            }
            
        } else if (j < crnt_ref_val->vb_cnt) { //if the current vbase is not the last vbase
            j++;                               //increment the vbase index.
        }
    }
    
    //the length test below is redundant with the length test performed during read mapping
    //but is kept here as a sanity check
    
    if (i != len) { //check that loop exited because the qscore for every base was tested
        pass = 0;   //if not, set pass to false
    }
    
    return pass;    //return the result of the test
}

/* count_matched_targets: count the number of targets to which at least 1 read mapped */
int count_matched_targets(target * trgts, target_params * trg_prms)
{
    int i = 0;               //general purpose index
    int matched_targets = 0; //number of targets with at least 1 mapped read
    
    for (i = 0; i < trg_prms->t_cnt; i++) {  //for every target
        if (!trgts[i].mul && trgts[i].cnt) { //if target is not redundant and has >=1 read mapped
            matched_targets++;               //increment matched targets count
        } 
    }
    return matched_targets;  //return matched targets count
}

/* crrct_testdata_nonsrc_mtch: decrement match counters when a mutant testdata read maps to a target
 from a different source sequence or maps to the correct target sequence by creating/extending a
 homopolymer at the 3' end of the target. */
void crrct_aberrant_testdata_mtch(uint32_t * match_cnt, opt_mx_trg * trg_vals, mapping_metrics * met, int channel, int chnl_mtch_typ, testdata_vars * testdata)
{
    //perform redundant check that testdata mode is ON
    if (!testdata->run) {
        printf("crrct_testdata_nonsrc_mtch: ERROR - non-source match correction can only be performed during testdata analysis. aborting...\n");
        abort();
    }

    (*match_cnt)--;                          //decrement target match count
    met->matches--;                          //decrement total match count
    met->chan_count[channel]--;              //decrement channel match count
    
    if (channel == BND) {                    //read is from the bound channel
        trg_vals->bnd--;                     //decrement target bound channel count
        
        if (chnl_mtch_typ == FULL) {         //read is a full bound channel barcode match
            met->full_match[BND]--;          //decrement full bound channel barcode matches
            
        } else if (chnl_mtch_typ == PART) {  //read is a part bound channel barcode match
            met->part_match[BND]--;          //decrement part bound channel barcode matches
        }
        
    } else if (channel == UNB) {             //read is from the unbound channel
        trg_vals->unb--;                     //decrement unbound channel count
        
        if (chnl_mtch_typ == FULL) {         //read is a full unbound channel barcode match
            met->full_match[UNB]--;          //decrement full unbound channel barcode match
            
        } else if (chnl_mtch_typ == PART) {  //read is a part unbound channel barcode match
            met->part_match[UNB]--;          //decrement part unbound channel barcode match
        }
    }
    
    return;
}

//
//  prcs_rds_DNA_qc.c
//  
//
//  Created by Eric Strobel on 9/8/26.
//

#include <stdio.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../../variant_maker/variant_maker_defs.h"
#include "../../../variant_maker/vmt_suffix.h"
#include "../../../variant_maker/make_barcodes.h"

#include "../../../TECdisplay_mapper/TECdisplay_mapper_defs.h"
#include "../../../TECdisplay_mapper/TECdisplay_mapper_structs.h"
#include "../../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"
#include "../../../TECdisplay_mapper/map_reads/map_expected/print_navigator_template.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "../UNV/bypass_fastp.h"
#include "../MUX/map_barcoded_targets.h"

#include "../../DNAQC_trgt_gen/mk_DNA_qc_trgts.h"
#include "../../DNAQC_trgt_gen/mk_DNA_qc_testdata.h"

#include "../MLT/prcs_MLT_cotrans.h"
#include "../MUX/get_brcd_str.h"

#include "../../DNAQC_trgt_gen/mk_DNA_qc_trgts.h"

#include "./print_DNA_qc_output.h"

#include "prcs_rds_DNA_qc.h"


struct testdata_DNA_qc_vars testdata_DNA_qc = {0}; //structure containing test data read analysis variables

/* prcs_rds_DNA_qc: manages processing of TECprobe-MUX data */
int prcs_rds_DNA_qc(TPROBE_names * nm, FILE * fp_MUXtrgs, int trgt_ftype, fastp_params fastp_prms, testdata_DNA_qc_vars * testdata_DNA_qc, int run_bypass_fastp)
{
    FILE *ifp[READ_MAX] = {NULL};  //pointers for input fastq files
    mapping_metrics met = {0};     //read processing metrics storage
    
    //TODO: this is probably not needed
    init_chnl_mtrcs_mem(&met, TPROBE_CHANNEL_MAX); //initialize channel tracking memory
    
    target * refs = {NULL};        //pointer for array of reference targets
    opt_ref * ref_val = {NULL};    //pointer for array of optional reference target structures
    
    compact_target * ctrg  = NULL; //pointer to compact target structures
    opt_BC * BC_val = NULL;        //pointer to optional target value structures
    association * assoc = NULL;    //pointer to association structure
    
    target_params trg_prms = {0};  //structure for storing target parameters
    TDSPLY_fasta wt = {0};         //storage for wt sequence information
    
    int ctrg_cnt = 0;              //number of compact targets stored
    int clcd_ctrg_cnt = 0;         //calculated number of barcode targets
    
    char line[MAX_LINE+1] = {0};   //array to store line
    
    //allocate memory for reference targets
    if ((refs = calloc(MAXREF, sizeof(*refs))) == NULL) {
        printf("prcs_rds_DNA_qc: error - reference target memory allocation failed\n");
        return 1;
    }
    
    if ((ref_val = calloc(MAXREF, sizeof(*ref_val))) == NULL) {
        printf("prcs_rds_DNA_qc: error - reference target value memory allocation failed\n");
        return 1;
    }
    
    if (trgt_ftype == FASTA_FILE) {
        
        //iterate through file to count number of targets
        while (get_line(line, fp_MUXtrgs)) { //until all lines have been read
            if (line[0] == '>') {            //if reading first line of fasta entry
                get_line(line, fp_MUXtrgs);  //get the second line of the fasta entry
                trg_prms.xpctd +=2;          //increment expected barcode count
                //NOTE: each barcode is assigned its own ctrg (i.e. each target seq gets two ctrgs)
            } else {
                printf("prcs_rds_DNA_qc: error - targets must be provided in fasta format. aborting...\n");
                abort();
            }
        }
        
        fclose(fp_MUXtrgs);                  //close targets file
        get_file(&(fp_MUXtrgs), nm->trgts);  //re-open targets file
        
    } else { //throw error if barcode target file type is not fasta
        printf("prcs_rds_DNA_qc: unrecognized barcoded target file type. aborting...\n");
        abort();
    }
    
    clcd_ctrg_cnt = trg_prms.xpctd * 129;  //caculate expected number of BC targs, including single subs and indels
    
    //allocate memory for compact targets
    if ((ctrg = calloc(clcd_ctrg_cnt, sizeof(*ctrg))) == NULL) {
        printf("prcs_rds_DNA_qc: error - compact target memory allocation failed. aborting...\n");
        abort();
    }
    
    //allocate memory for barcode target optional values
    if ((BC_val = calloc(clcd_ctrg_cnt, sizeof(*BC_val))) == NULL) {
        printf("prcs_rds_DNA_qc: error - optional target value memory allocation failed. aborting...\n");
        abort();
    }
    
    //allocate memory for association structures
    if ((assoc = calloc(trg_prms.xpctd, sizeof(*assoc))) == NULL) {
        printf("prcs_rds_DNA_qc: error - association memory allocation failed. aborting...\n");
        abort();
    }
        
    //generate barcode targets
    ctrg_cnt = mk_DNA_qc_trgts(refs, ref_val, ctrg, BC_val, assoc, fp_MUXtrgs, trgt_ftype, &trg_prms, clcd_ctrg_cnt, &wt);
    met.srcTrgs = trg_prms.t_cnt; //record number of source barcodes
    met.targets = ctrg_cnt;       //record number of targets generated
    
    printf("src: %d\ntrg: %d\n", trg_prms.t_cnt, ctrg_cnt);
    
    /********* hash table initialization and construction **********/
    compact_h_node **htbl_MUX = NULL;                  //hash table root
    compact_h_node_bank hn_MUX_bank = {NULL, NULL, 0}; //bank for hash table nodes
    
    hn_MUX_bank.count = 0; //initialize hash table node count to zero
    
    //allocate TABLE_SIZE hash table node pointers
    if ((htbl_MUX = calloc(TABLE_SIZE, sizeof(*htbl_MUX))) == NULL) {
        printf("prcs_rds_DNA_qc: error - hash table memory allocation failed\n");
        abort();
    }
    
    //allocate BLOCK_SIZE hash table nodes
    if ((hn_MUX_bank.chn = calloc(BLOCK_SIZE, sizeof(*(hn_MUX_bank.chn)))) == NULL) {
        printf("prcs_rds_DNA_qc: error - hash table node memory allocation failed\n");
        abort();
    }
    
    compact_h_node_bank *crrnt_hn_MUX_bank = &hn_MUX_bank;                     //pointer for handling hash node bank
    mk_htbl_MUX(htbl_MUX, crrnt_hn_MUX_bank, ctrg, ctrg_cnt, &trg_prms, &met); //generate barcode target hash table
    /****** end of hash table initialization and construction ******/
    
    /*************** testdata generation *****************/
    if (testdata_DNA_qc->run) {
        mk_DNA_qc_testdata(nm, ctrg, ctrg_cnt);
    }
    /*********** end of test data generation *************/
    
    char fp_out_dir[10] = {"fastp_out"};
    mk_out_dir(fp_out_dir); //make directory for output files
    
    //perform or bypass fastp processing
    if (testdata_DNA_qc->run && run_bypass_fastp) {
        bypass_fastp(nm->file[READ1], nm->file[READ2], &ifp[0], fp_out_dir);
    } else {
        call_fastp_TPROBE(nm->file[READ1], nm->file[READ2], &ifp[0], fastp_prms, fp_out_dir); //fastp pre-processing
    }
    
    int pairing[2]; //variable for tracking number of concordant (index 0) and discordant (index 1) reads
        
    assess_brcd_concordance(&ifp[0], htbl_MUX, nm, ctrg, trg_prms.t_cnt, ctrg_cnt, &pairing[0], &met, fastp_prms.mode);
    
    if (testdata_DNA_qc->run) {                        //if running testdata mode
        assess_testdata_mapping(ctrg, trg_prms.t_cnt); //assess whether testdata mapped as expected
    }
    
    print_DNA_qc_output(NULL, ctrg, trg_prms.t_cnt, &pairing[0]); //print output
    
    return 1;
}

/* assess_brcd_concordance: determine whether barcodes present in a sequencing read are concordant */
void assess_brcd_concordance(FILE **ifp, compact_h_node **htbl_MUX, TPROBE_names * nm, compact_target * ctrg, int brcd_cnt, int ctrg_cnt, int * pairing, mapping_metrics * met, int mode)
{
    extern int debug;                                   //flag to turn on debug mode
    extern struct testdata_DNA_qc_vars testdata_DNA_qc; //structure containing test data read analysis variables
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    //general variables
    int got_line[READ_MAX] = {0};   //flag to indicate success of the get_line function
    int proceed = 1;                //flag to indicate that read processing should proceed
    
    //read variables
    char read1[FQ_LINES][MAX_LINE] = {{0}};  //array for storing all four lines of one read 1 fastq entry
    char read2[FQ_LINES][MAX_LINE] = {{0}};  //array for storing all four lines of one read 2 fastq entry
    char ipt_ID[READ_MAX][MAX_LINE] = {{0}}; //array for storing input lines when debug mode is on
    char ipt[READ_MAX][MAX_LINE] = {{0}};    //array for storing input reads for error checking
    
    //barcode variables
    char brcd_str1[MAX_BARCODE_LEN+1] = {0};    //barcode string 1
    char brcd_str2[MAX_BARCODE_LEN+1] = {0};    //barcode string 2
    char rc_brcd_str2[MAX_BARCODE_LEN+1] = {0}; //reverse complement of barcode string 1
    
    int r1_printed = -1; //flag to indicate that read 1 printing was successful
    int r2_printed = -1; //flag to indicate that read 2 printing was successful
    
    compact_target * crnt_mpd_trg1 = NULL; //pointer to mapped target
    compact_target * crnt_mpd_trg2 = NULL; //pointer to mapped target
    
    compact_target * crnt_ref_trg1 = NULL; //pointer to reference target of mapped target
    compact_target * crnt_ref_trg2 = NULL; //pointer to reference target of mapped target
    
    opt_BC * BC_val1 = NULL; //pointer to barcode target optional values
    opt_BC * BC_val2 = NULL; //pointer to barcode target optional values
    
    association * bc1_assoc = NULL; //pointer to association structure
    association * bc2_assoc = NULL; //pointer to association structure
    
    int bc1_cncrdnt = 0; //flag that barcode 1 is concordant with barcode 2
    int bc2_cncrdnt = 0; //flag that barcode 2 is concordant with barcode 1
        
    for (i = 0; proceed; i++) {
        
        brcd_str1[0] = brcd_str2[0] = '\0';   //initialize barcode strings to 0
        crnt_mpd_trg1 = crnt_mpd_trg2 = NULL; //initialize reference and mapped target pointers to NULL
        bc1_cncrdnt = bc2_cncrdnt = 0;        //initialize concordance flags to 0
        
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
                strcpy(ipt_ID[READ2], read2[LINE1]); //store read 2 ID for post-processing comparison
            }
            
            if (!(got_line[READ1] && got_line[READ2])) { //test success of get_line for read 1 and read 2
                //test failed. this should only occur at end of file and get_line will
                //throw an error if an unexpected input line is encountered
                
                //test that the end of read 1 and 2 fastq files were reached at the same time
                if (got_line[READ1] == got_line[READ2]) { //reached EOF for both files
                    proceed = 0; //exit read processing loop
                } else {
                    printf("assess_brcd_concordance: error R1out.fq and R2out.fq (fastp output) do not contain the same number of lines. aborting...\n");
                    abort();
                }
            }
        }
        
        verify_read(read1); //verify read 1 integrity
        verify_read(read2); //verify read 2 integrity
        
        //identify barcodes and assess concordance
        if (proceed) {
            met->reads_processed++;                //track total number of reads processed
            
            get_brcd_str(brcd_str2, brcd_str1, 1, 1, &read1[LINE1][0]); //get barcode from read1 id
            reverse_complement(rc_brcd_str2, brcd_str2, REVCOMP);       //revcomp barcode 2 string
                        
            crnt_ref_trg1 = map_brcd(brcd_str1, htbl_MUX, &crnt_mpd_trg1, met);    //map barcode1 using hash table
            crnt_ref_trg2 = map_brcd(rc_brcd_str2, htbl_MUX, &crnt_mpd_trg2, met); //map barcode2 using hash table
            
            if (crnt_ref_trg1 != NULL && crnt_ref_trg2 != NULL) { //if barcodes mapped
                
                if (testdata_DNA_qc.run) { //if mapping test data, check that read mapped to correct target
                    //compare_testdata_barcode_id(crnt_mpd_trg, &read2[LINE1][0], &read2[LINE2][0]);
                    //TODO: determine whether to include this test
                }
                
                BC_val1 = (opt_BC *)crnt_ref_trg1->opt; //set pointer to ref barcode optional values
                BC_val2 = (opt_BC *)crnt_ref_trg2->opt; //set pointer to ref barcode optional values
                
                bc1_assoc = (association *)(crnt_ref_trg1->utl); //set barcode 1 association structure
                bc2_assoc = (association *)(crnt_ref_trg2->utl); //set barcode 2 association structure
                
                //test barcode 1 and barcode 2 concordance separately
                bc1_cncrdnt = ((uint64_t)(BC_val1->lnk) == (uint64_t)(crnt_ref_trg2)) ? 1 : 0;
                bc2_cncrdnt = ((uint64_t)(BC_val2->lnk) == (uint64_t)(crnt_ref_trg1)) ? 1 : 0;
                
                if (bc1_cncrdnt && bc2_cncrdnt) { //barcodes 1 and 2 are concordant
                    pairing[CONCORDANT]++; //increment total number of concordant reads
                    bc1_assoc->cncrdnt++;  //increment number of barcode 1 concordant reads
                    bc2_assoc->cncrdnt++;  //increment number of barcode 2 concordant reads
                    
                } else if (!bc1_cncrdnt && !bc2_cncrdnt) { //barcodes 1 and 2 are discordant
                    pairing[DISCORDANT]++; //increment total number of discordant reads
                    bc1_assoc->dscrdnt++;  //increment number of barcode 1 discordant reads
                    bc2_assoc->dscrdnt++;  //increment number of barcode 2 discordant reads
                    
                } else { //if concordance is inconsistent, throw error and abort
                    printf("assess_brcd_concordance: error - inconsistent concordance. aborting...\n");
                    abort();
                }
                                
                // *** verify read sequences ***

                //this is a sanity check, there is no reason
                //that the read sequences should ever change
                            
                if (!strcmp(ipt[READ1], read1[LINE2])) {
                    met->read_matches[READ1]++;
                } else {
                    printf("assess_brcd_concordance: CRITICAL - sequence integrity failure in read 1. this should NEVER happen! please contact estrobel@buffalo.edu regarding this error.\naborting...\n");
                    abort();
                }
                
                if (!strcmp(ipt[READ2], read2[LINE2])) {
                    met->read_matches[READ2]++;
                } else {
                    printf("assess_brcd_concordance: CRITICAL - sequence integrity failure in read 2. this should NEVER happen! please contact estrobel@buffalo.edu regarding this error.\naborting...\n");
                    abort();
                }
                // *** end verify read sequences ***
                
                //track ongoing processing
                if (i+1 % 1000000 == 0) {
                    printf(">");
                }
            }
        }
    }
    
    /*
    for (i = 0; i < brcd_cnt; i++) {
        bc1_assoc = (association *)(ctrg[i].utl);
        printf("%d\t%d\t%d\t%d\t%d\t%d\n", i, bc1_assoc->cncrdnt, bc1_assoc->dscrdnt, bc1_assoc->src_cncrdnt, bc1_assoc->src_dscrdnt, bc1_assoc->rnd_dscrdnt);
    }
    */
    
    return;
}

/* assess_testdata_mapping: assess whether testdata mapping matches expectations */
void assess_testdata_mapping(compact_target * ctrg, int brcd_cnt)
{
    int i = 0; //general purpose index
    
    opt_BC * BC_val = NULL; //pointer to barcode values structure
    
    compact_target * lnkd_ctrg = NULL; //pointer to linked compact target
    
    association * crnt_assoc = NULL; //pointer to current target association structure
    association * lnkd_assoc = NULL; //pointer to linked target association structure
    
    int fail = 0; //flag that testdata analysis failed
    
    for (i = 0; i < brcd_cnt; i++) {      //for every barcode target
        
        BC_val = (opt_BC *)(ctrg[i].opt); //set pointer to current target barcode values
        crnt_assoc = ctrg[i].utl;         //set pointer to current target association structure
        
        lnkd_ctrg = BC_val->lnk;          //set pointer to linked target barcode values
        lnkd_assoc = lnkd_ctrg->utl;      //set pointer to linked target association structure
        
        //if the number of concordant mapped reads does not equal the number of concordant testdata
        //reads generated with the current target and its linked target as the source barcode,
        //testdata analysis fails.
        if (crnt_assoc->cncrdnt != (crnt_assoc->src_cncrdnt + lnkd_assoc->src_cncrdnt)) {
            printf("assess_testdata_mapping: error - unexpected number of concordant barcodes\n");
            printf("expected: %d\n", crnt_assoc->src_cncrdnt + lnkd_assoc->src_cncrdnt);
            printf("  mapped: %d\n", crnt_assoc->cncrdnt);
            printf("current barcode: %s\tsrc_cncrdnt: %d\n", ctrg[i].cid, crnt_assoc->src_cncrdnt);
            printf(" linked barcode: %s\tsrc_cncrdnt: %d\n\n", lnkd_ctrg->cid, lnkd_assoc->src_cncrdnt);
            fail = 1;
        }
        
        //if the number of discordant mapped reads does not equal the number of discordant testdata
        //reads generated with the current target as the source plus the numer of discordant
        //testdata reads in which the current target barcode was randomly selected as the
        //discordant barcode, testdata analysis fails
        if (crnt_assoc->dscrdnt != (crnt_assoc->src_dscrdnt + crnt_assoc->rnd_dscrdnt)) {
            printf("assess_testdata_mapping: error - unexpected number of discordant barcodes\n");
            printf("expected: %d\n", crnt_assoc->src_dscrdnt + crnt_assoc->rnd_dscrdnt);
            printf("  mapped: %d\n", crnt_assoc->dscrdnt);
            printf(" barcode: %s\tsrc_dscrdnt: %d\trnd_dscrdnt: %d\n\n", ctrg[i].cid, crnt_assoc->src_dscrdnt, crnt_assoc->rnd_dscrdnt);
            fail = 1;
        }
    }
    
    //print pass/fail message
    if (fail) {
        printf("\ntestdata analysis failed. aborting...\n\n");
        abort();
    } else {
        printf("\ntestdata_analysis passed\n\n");
    }
    
    return;
}

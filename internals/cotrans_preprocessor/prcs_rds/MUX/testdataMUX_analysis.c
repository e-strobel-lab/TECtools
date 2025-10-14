//
//  testdataMUX_analysis.c
//  
//
//  Created by Eric Strobel on 7/10/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/gen_utils.h"
#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "../../MUX_trgt_gen/mk_MUX_trgts.h"

#include "../../../variant_maker/make_barcodes.h"

#include "testdataMUX_analysis.h"

struct testdata_MUX_vars testdata_MUX = {0}; //structure containing test data read analysis variables

/* get_testdata_barcode_id: get barcode id from test data read line 1 */
uint64_t get_testdata_barcode_id(char * id)
{
    extern struct testdata_MUX_vars testdata_MUX;
    
    char * p_bid = NULL;            //pointer to barcode id in read id line
    char bid_str[MAX_LINE+1] = {0}; //array for storing barcode id string
    
    //set pointer to string barcode=' in testdata read id
    if ((p_bid = strstr(id, "barcode=")) == NULL) {
        printf("TEST DATA ANALYSIS: error - read line 1 does not contain expected barcode information. aborting...\n");
        abort();
    }
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    //store barcode id string
    for (i = 0, j = strlen("barcode="); isdigit(p_bid[j]) && i < MAX_UINT64_LEN; i++, j++) {
        bid_str[i] = p_bid[j];
    }
    bid_str[i] = '\0';
    
    //check that barcode id was formatted correctly in read id
    if (p_bid[j] != '_' || i < MAX_UINT64_LEN) {
        printf("get_testdata_barcode_id: unexpected format for barcode id in test data read. aborting...\n");
        abort();
    }
    
    //return barcode id value
    return strtoull(bid_str, NULL, 10);
}

/* compare_testdata_ barcode id: check whether test data read mapped to the expected barcode */
void compare_testdata_barcode_id(compact_target * ctrg, char * rd_id, char * sq)
{
    extern struct testdata_MUX_vars testdata_MUX;
    
    uint64_t bid = get_testdata_barcode_id(rd_id); //get barcode id from test data read line 1
    
    //test whether barcode id of target that read mapped to matches the expected barcode id
    if (bid != ctrg->bid) {
        printf("ERROR: test data read did not map to the correct barcode or barcode mutant:\n");
        printf("read id: %s\n", rd_id);
        printf("brcd id: %llu\n", (long long unsigned int)bid);
        printf("ctrg id: %llu\n", (long long unsigned int)ctrg->bid);
        printf("read2:    %s\n", sq);
    } else {
        testdata_MUX.bid_match++;
    }
}


/* print_MUX_testdata_analysis: report outcome of test data read mapping */
void print_MUX_testdata_analysis(mapping_metrics * met, compact_target * ctrg)
{
    extern struct testdata_MUX_vars testdata_MUX;
    
    int i = 0;                    //general purpose index
    char out_str[MAX_LINE] = {0}; //string for storing output line
    
    char pass[7] = "[PASS]";  //output string for PASS outcomes
    char fail[7] = "[FAIL]";  //output string for FAIL outcomes
    char * pf_ptr = NULL;     //pointer to PASS/FAIL strings
    
    //open testdata analysis output file
    FILE * testdata_fp = NULL;
    if ((testdata_fp = fopen("test_data_analysis.txt", "w")) == NULL) {
        printf("print_MUX_testdata_analysis: ERROR - could not generate test data analysis file. aborting ...\n");
        abort();
    }
    
    printf("\n");
    sprintf(out_str, "***** TEST DATA ANALYSIS SUMMARY *****\n");
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report that test data reads contain TECprobe-MUX barcodes
    sprintf(out_str, ">    mode - sequences contain TECprobe-MUX barcodes\n");
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report number of processed test data reads
    sprintf(out_str, ">   depth - %d sequences were analyzed\n", met->reads_processed);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report barcode mapping efficiency
    pf_ptr = (met->mapped == met->reads_processed) ? &pass[0] : &fail[0];
    sprintf(out_str, "> barcode - %s %10.6f%% (expected 100%%) of sequences mapped to a barcode target\n",
            pf_ptr, ((float)(met->mapped)/(float)(met->reads_processed))*100);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report barcode mapping accuracy
    pf_ptr = (testdata_MUX.bid_match == met->mapped) ? &pass[0] : &fail[0];
    sprintf(out_str, "> barcode - %s %10.6f%% (expected 100%%) of mapped sequences mapped to the expected barcode target\n",
            pf_ptr, ((float)(testdata_MUX.bid_match)/(float)(met->mapped))*100);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report channel mapping
    double test_val = 0; //used to store calculated values
    
    double expctd_chan[3] = {0.4, 0.4, 0.2};
    char chan_nm[3][32] = {"untreated", "modified", "unmappable"};
    
    for (i = 0; i < TPROBE_CHANNEL_MAX; i++) {
        //calculate fraction of reads in channel i
        test_val = (double)(met->chan_count[i])/(double)(met->reads_processed);
        
        //check if fraction of reads in channel i matches expected value
        pf_ptr = (compare_float(test_val, expctd_chan[i], 0.000001)) ? &pass[0] : &fail[0];
        
        //print output
        sprintf(out_str, "> channel - %s %10.6f%% (expected  %.0f%%) of sequences contained a(n) %s barcode\n",
                pf_ptr, test_val*100, expctd_chan[i]*100, chan_nm[i]);
        printf2_scrn_n_fl(testdata_fp, out_str);
    }
    
    //report percentage of full and partial barcode matches for untreated and modified channels
    for (i = 0; i < 2 ; i++) {
        //report percentage of reads with a full barcode match
        test_val = (double)(met->full_match[i])/(double)(met->chan_count[i]);
        pf_ptr = (compare_float(test_val, 0.5, 0.000001)) ? &pass[0] : &fail[0];
        sprintf(out_str, "> channel - %s %10.6f%% (expected  50%%) of %s sequences contained a full barcode match\n",
                pf_ptr, test_val*100, chan_nm[i]);
        printf2_scrn_n_fl(testdata_fp, out_str);
        
        //report percentage of reads with a partial barcode match
        test_val = (double)(met->part_match[i])/(double)(met->chan_count[i]);
        pf_ptr = (compare_float(test_val, 0.5, 0.000001)) ? &pass[0] : &fail[0];
        sprintf(out_str, "> channel - %s %10.6f%% (expected  50%%) of %s sequences contained a partial barcode match\n",
                pf_ptr, test_val*100, chan_nm[i]);
        printf2_scrn_n_fl(testdata_fp, out_str);
    }
    
    if ((fclose(testdata_fp)) == EOF) {
        printf("print_MUX_testdata_analysis: error - error occurred when closing test data analysis file. Aborting program...\n");
        abort();
    }
}

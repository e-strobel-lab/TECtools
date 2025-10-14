//
//  testdata3pEnd_analysis.c
//  
//
//  Created by Eric Strobel on 3/15/22.
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
#include "../../../seq_utils/mapping_metrics.h"

#include "testdata3pEnd_analysis.h"


struct testdata_3pEnd_vars testdata_3pEnd = {0}; //structure containing test data read analysis variables

/* check_testdata_ipt: verify testdata analysis input files */
int check_testdata_ipt(TPROBE_names * nm)
{
    //when using test data, check whether the supplied reads are actually test data fastq files and
    //whether the test data reads contain native or randomized 3' ends
    //test data filenames must always contain the string "testdata<MODE>_R*.fq", where <MODE> is 'NATIVE'
    //or 'RANDOM' to indicate the read generation mode and '*' is 1 or 2 to indicate read 1 and 2 files
    if (testdata_3pEnd.run) {
        if (strstr(nm->file[READ1], "testdataNATIVE_R1.fq") != NULL &&
            strstr(nm->file[READ2], "testdataNATIVE_R2.fq") != NULL) {
            testdata_3pEnd.mode = TESTDATA_NAT;
        } else if (strstr(nm->file[READ1], "testdataRANDOM_R1.fq") != NULL &&
                   strstr(nm->file[READ2], "testdataRANDOM_R2.fq") != NULL) {
            testdata_3pEnd.mode = TESTDATA_RND;
        } else {
            printf("TEST DATA ANALYSIS: error - test data filename does not have expected NATIVE or RANDOM identifier string. aborting...\n");
            abort();
        }
    }
    
    //next, check that read 1, read 2, ends files were generated using the same end target length
    //in reads filenames, a two digit end target length value will precede the string "ntLong_testdata"
    //in the ends filename, a two digit end target length value will precede the string "ntLong.txt"
    char * trg_len_ptr = NULL;			//pointer to end target length value
    char tmp_trg_len[3] = {0};			//array for storing end target length characters
    int trg_len_rds[READ_MAX] = {0};	//end target length from reads filenames
    int trg_len_ends = 0;				//end target length from the ends filename
    
    int i = 0;
    
    //get end target length values from reads 1 and 2 filenames
    for (i = 0; i < READ_MAX; i++) {
        
        //set pointer to start of string that follows two-digit end target length value
        if ((trg_len_ptr = strstr(nm->file[i], "ntLong_testdata")) != NULL) {
            
            //copy end target length from filenames to tmp_trg_len array
            tmp_trg_len[0] = trg_len_ptr[-2];
            tmp_trg_len[1] = trg_len_ptr[-1];
            tmp_trg_len[2] = '\0';
            
            //check that end target length characters are digits
            if (!isdigit(tmp_trg_len[0]) && !isdigit(tmp_trg_len[1])) {
                printf("TEST DATA ANALYSIS: error - unexpected target length format in testdata filename. aborting...\n");
                abort();
            }
            
            trg_len_rds[i] = atoi(tmp_trg_len); //store end target length value
        } else {
            printf("TEST DATA ANALYSIS: error - test data read %d filename does not contain expected target length string. aborting...\n", i+1);
            abort();
        }
    }
    
    //get end target length values from ends targets filename
    //set pointer to start of string that follows two-digit end target length value
    if ((trg_len_ptr = strstr(nm->trgts, "ntLong.txt")) != NULL) {
        
        //copy end target length from filenames to tmp_trg_len array
        tmp_trg_len[0] = trg_len_ptr[-2];
        tmp_trg_len[1] = trg_len_ptr[-1];
        tmp_trg_len[2] = '\0';
        
        //check that end target length characters are digits
        if (!isdigit(tmp_trg_len[0]) && !isdigit(tmp_trg_len[1])) {
            printf("TEST DATA ANALYSIS: error - unexpected target length format in end targets filename. aborting...\n");
            abort();
        }
        
        trg_len_ends = atoi(tmp_trg_len); //store end target length value
    } else {
        printf("TEST DATA ANALYSIS: error - 3' end targets filename does not contain expected target length string. aborting...\n");
        abort();
    }
    
    //check that both read 1 and read 2 filenames contain an end target length value
    //matching that of the ends targets filename
    if ((trg_len_rds[READ1] != trg_len_ends) || (trg_len_rds[READ2] != trg_len_ends)) {
        printf("TEST DATA ANALYSIS: error - target length strings in input sequence and 3' end targets filenames do not match. aborting...\n");
        abort();
    }
    
    return 1;
}


/* get_testdata_3pEND: get 3' end info from testdata read id */
void get_testdata_3pEND(char * id)
{
    extern struct testdata_3pEnd_vars testdata_3pEnd; //structure containing test data read analysis variables
    
    testdata_3pEnd.end3p = 0;	//zero 3' end value
    
    char * ptr_3pEnd = NULL;	//pointer to string '3pEnd=' in test data read id
    char end_cha[4] = {0};		//3' end value string
    
    //set pointer to string 3pEnd=' in testdata read id
    if ((ptr_3pEnd = strstr(id, "3pEnd=")) == NULL) {
        printf("TEST DATA ANALYSIS: error - read line 1 does not contain expected 3' end information. aborting...\n");
        abort();
    }
    
    //set end_char to 3' end value in testdata read id
    //test data 3' end value will always be a three digit number
    int i = 0;
    int j = 0;
    
    //j is initialized to six to offset ptr_3pEnd to the
    //character immediately after the string "3pEND="
    for (i = 0, j = 6; i < 3; i++, j++) {
        if (isdigit(ptr_3pEnd[j])) {
            end_cha[i] = ptr_3pEnd[j];
        } else {
        	printf("TEST DATA ANALYSIS: error - expected 3' end value in read line 1 contains a non-digit character. aborting...\n");
            abort();
        }
    }
    end_cha[i] = '\0';
    
    testdata_3pEnd.end3p = atoi(end_cha); //store 3' end value as int in testdata structure
}

/* compare_testdata_3pEnd: compare expected 3' end value, which is contained in line 1 of testdata
 fastq entries, with thes mapped 3' end value */

void compare_testdata_3pEnd(int end, char id[MAX_LINE], char sq[MAX_LINE])
{
    extern struct testdata_3pEnd_vars testdata_3pEnd; //structure containing test data read analysis variables
    
    int distance = 0;
    
    if (end) { //end mapped
        get_testdata_3pEND(&id[0]);			//get expected 3' end from read id
        if (testdata_3pEnd.end3p == end) {	//mapped 3' end matches expected value
            testdata_3pEnd.perfect[end]++; 	//count exact matches to expected 3' end
        }
        
        distance = abs(testdata_3pEnd.end3p-end); //calculate distance between expected and mapped end
        testdata_3pEnd.distance[distance]++;	  //track distance of mapped 3' end from expected value
        
        //track maximum observed distance between expected and mapped 3' end values
        if (distance > testdata_3pEnd.max_distance) {
            testdata_3pEnd.max_distance = distance;
        }
        
    } else { //end did not map
        //print the id and read 2 sequence of the unmapped read
        //this follows the TEST DATA ANALYSIS output from the map_3pEnd function,
        //which reports the seed sequence of the read that did not map. this
        //information is crucial because every test data read should be mappable
        printf("read 2 id: %s\nread 2 sq: %s\n", id, sq);
    }
}

/* print_3pEnd_testdata_analysis: assess testdata analysis
 outcomes and print results to a file and to the screen */
void print_3pEnd_testdata_analysis(mapping_metrics *met, target3p_params trg_prms, target * trgts)
{
    extern struct testdata_3pEnd_vars testdata_3pEnd; //structure containing test data read analysis variables
    
    int i = 0;
    char out_str[8192] = {0}; //string for storing output line
    
    char pass[7] = "[PASS]";  //output string for PASS outcomes
    char fail[7] = "[FAIL]";  //output string for FAIL outcomes
    char * pf_ptr = NULL;	  //pointer to PASS/FAIL strings
    
    //open testdata analysis output file
    FILE * testdata_fp = NULL;
    if ((testdata_fp = fopen("test_data_analysis.txt", "w")) == NULL) {
        printf("print_3pEnd_testdata_analysis: ERROR - could not generate test data analysis file. aborting ...\n");
        abort();
    }
    
    printf("\n");
    sprintf(out_str, "***** TEST DATA ANALYSIS SUMMARY *****\n");
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report whether input test data reads are expected to map to native or random targets
    switch (testdata_3pEnd.mode) {
        case TESTDATA_NAT: sprintf(out_str, ">    mode - sequences contain NATIVE 3' ends\n"); break;
        case TESTDATA_RND: sprintf(out_str, ">    mode - 3' ends contain RANDOM 1nt substitutions and indels\n"); break;
        default:
            break;
    }
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report number of processed test data reads
    sprintf(out_str, ">   depth - %d sequences were analyzed\n", met->reads_processed);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report 3' end mapping efficiency
    pf_ptr = (met->mapped == met->reads_processed) ? &pass[0] : &fail[0];
    sprintf(out_str, ">  3' end - %s %10.6f%% (expected 100%%) of sequences mapped to a 3' end target\n",
            pf_ptr, ((float)(met->mapped)/(float)(met->reads_processed))*100);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report native target mapping efficiency when using native test data reads
    if (testdata_3pEnd.mode == TESTDATA_NAT) {
        
        //report whether all test data reads mapped to native 3' end targets
        pf_ptr = (met->nat_cnt == met->reads_processed) ? &pass[0] : &fail[0];
        sprintf(out_str, ">  3' end - %s %10.6f%% (expected 100%%) of sequences mapped to a native 3' end target\n",
                pf_ptr, ((float)(met->nat_cnt)/(float)(met->reads_processed))*100);
        printf2_scrn_n_fl(testdata_fp, out_str);
        
        //report whether all test data reads mapped with a 0nt distance from the expected target
        pf_ptr = (testdata_3pEnd.max_distance == 0) ? &pass[0] : &fail[0];
        sprintf(out_str, ">  3' end - %s the max distance between mapped and expected transcript lengths is %d (expected 0)\n", pf_ptr, testdata_3pEnd.max_distance);
        printf2_scrn_n_fl(testdata_fp, out_str);
        
        if (pf_ptr[1] == 'F') { //check if max distance test failed, only print this to screen
            printf(">  3' end - [NOTE] see test_data_analysis.txt for a 3' end mapping distance table\n");
        }
    }
    
    //report 3' end target coverage when using 1 nt variant test reads
    if (testdata_3pEnd.mode == TESTDATA_RND) {
        
        //report 3' end target coverage
        for (i = 0, pf_ptr = &pass[0]; i < trg_prms.cnt; i++) {
            //check whether reads mapped to each target
            //targets with mul=FALSE were the first instance of non-unique targets
            //and were therefore entered into the hash table and used to map reads
            if (!trgts[i].mul && !trgts[i].cnt) {
                sprintf(out_str, ">  3' end - [WARNING] zero sequences mapped to target %s\n", trgts[i].id);
                printf2_scrn_n_fl(testdata_fp, out_str);
                pf_ptr = &fail[0];
            }
        }
        
        
        if (pf_ptr[1] == 'P') {
            sprintf(out_str, ">  3' end - %s every 3' end target was used when mapping sequences\n", pf_ptr);
        } else {
            sprintf(out_str, ">  3' end - %s incomplete 3' end target coverage. try using a larger test data set\n", pf_ptr);
        }
        printf2_scrn_n_fl(testdata_fp, out_str);
        
        //report whether all test data reads mapped with a <=2 nt distance from the expected target
        pf_ptr = (testdata_3pEnd.max_distance <= 2) ? &pass[0] : &fail[0];
        sprintf(out_str, ">  3' end - %s the max distance between mapped and expected transcript lengths is %d (expected <=2)\n", pf_ptr, testdata_3pEnd.max_distance);
        printf2_scrn_n_fl(testdata_fp, out_str);
        
        //only print to screen
        printf(">  3' end - [NOTE] see test_data_analysis.txt for a 3' end mapping distance table\n");
        
    }
    
    double test_val = 0; //used to store calculated values
    
    double expctd_chan[3] = {0.45, 0.45, 0.1};
    char chan_nm[3][32] = {"untreated", "modified", "unmappable"};
    
    for (i = 0; i < TPROBE_CHANNEL_MAX; i++) {
        //calculate fraction of reads in channel i
        test_val = (double)(met->chan_count[i])/(double)(met->reads_processed);
        
        //check if fraction of reads in channel i matches expected value
        pf_ptr = (compare_float(test_val, expctd_chan[i], 0.000001)) ? &pass[0] : &fail[0];
        
        //print output
        sprintf(out_str, "> channel - %s %10.6f%% (expected  %f%%) of sequences contained an %s barcode\n",
                pf_ptr, test_val*100, expctd_chan[i]*100, chan_nm[i]);
        printf2_scrn_n_fl(testdata_fp, out_str);
    }


    //report percentage of full and partial barcode matches for untreated and modified channels
    for (i = 0; i < 2 ; i++) {
        //report percentage of reads with a full barcode match
        test_val = (double)(met->full_match[i])/(double)(met->chan_count[i]);
        pf_ptr = (compare_float(test_val, 0.9, 0.000001)) ? &pass[0] : &fail[0];
        sprintf(out_str, "> channel - %s %10.6f%% (expected  90%%) of %s sequences contained a full barcode match\n",
                pf_ptr, test_val*100, chan_nm[i]);
        printf2_scrn_n_fl(testdata_fp, out_str);
        
        //report percentage of reads with a partial barcode match
        test_val = (double)(met->part_match[i])/(double)(met->chan_count[i]);
        pf_ptr = (compare_float(test_val, 0.1, 0.000001)) ? &pass[0] : &fail[0];
        sprintf(out_str, "> channel - %s %10.6f%% (expected  10%%) of %s sequences contained a partial barcode match\n",
                pf_ptr, test_val*100, chan_nm[i]);
        printf2_scrn_n_fl(testdata_fp, out_str);
    }
    
    
    //print the distribution of mapped vs. expected end distances
    //(only in test_data_analysis.txt file
    fprintf(testdata_fp, "TEST DATA ANALYSIS: distance table\n");
    for (i = 0; i <= testdata_3pEnd.max_distance && i <= trg_prms.max; i++) {
        fprintf(testdata_fp, "%3d:\t%7d (%6.2f%%)\n", i, testdata_3pEnd.distance[i],
                ((float)(testdata_3pEnd.distance[i])/(float)(met->reads_processed))*100);
    }
    fprintf(testdata_fp, "distance >%d was not observed\n", testdata_3pEnd.max_distance);
    
    if ((fclose(testdata_fp)) == EOF) {
        printf("print_3pEnd_testdata_analysis: error - error occurred when closing test data analysis file. Aborting program...\n");
        abort();
    }
}

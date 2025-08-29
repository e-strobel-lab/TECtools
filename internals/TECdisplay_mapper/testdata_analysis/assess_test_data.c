//
//  assess_test_data.c
//  
//
//  Created by Eric Strobel on 5/3/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../TECdisplay_mapper_structs.h"
#include "../TECdisplay_mapper_defs.h"

#include "../../utils/gen_utils.h"

#include "assess_test_data.h"

/* parse_testdata_id: parse test data id line and store attributes*/
void parse_testdata_id(testdata_vars * testdata, char **td_trg_id, int * crnt_mut_cd, char * id_line)
{
    char * eq_loc = NULL;         //pointer to equals sign location in test data id line
    char * td_mut_cd = NULL;      //pointer to mutation code in test data id line
    
    //set start of target id
    if ((eq_loc = strstr(id_line, "=")) == NULL) { //locate equals sign in test data read id
        printf("map_expected_reads: error - unrecognized testdata read id format. aborting\n");
        abort();
    }
    
    if (isdigit(eq_loc[1]) || eq_loc[1] == '-') {  //check that equals sign is followed by a digit or '-'
        *td_trg_id = &eq_loc[1]; //set start of target id
    } else {
        printf("map_expected_reads: error - unrecognized testdata read id format. aborting\n");
        abort();
    }
    
    //set mutation code
    if ((td_mut_cd = strstr(id_line, "NAT")) != NULL) {        //test data read contains native sequence
        testdata->nat_rds++;                                   //increment native test data read count
        *crnt_mut_cd = NAT;                                    //set current mutation code to native
        
    } else if ((td_mut_cd = strstr(id_line, "SUB")) != NULL) { //test data read is a substitution mutant
        testdata->mut_rds++;                                   //increment mutant test data read count
        *crnt_mut_cd = SUB;                                    //set current mutation code to substitution
        
    } else if ((td_mut_cd = strstr(id_line, "INS")) != NULL) { //test data read is an insertion mutant
        testdata->mut_rds++;                                   //increment mutant test data read count
        *crnt_mut_cd = INS;                                    //set current mutation code to insertion
        
    } else if ((td_mut_cd = strstr(id_line, "DEL")) != NULL) { //test data read is a deletion mutant
        testdata->mut_rds++;                                   //increment mutant test data read count
        *crnt_mut_cd = DEL;                                    //set current mutation code to deletion
        
    } else {  //unrecognized mutation code
        printf("map_expected_reads: error - could not find mutation code in read id. aborting\n");
        abort();
    }
    
    //set end of target id
    if (td_mut_cd[-1] == '_') { //mutation code should be preceded by a '_'
        td_mut_cd[-1] = '\0';   //terminate target id string
    } else {
        printf("map_expected_reads: error - unrecognized testdata read id format. aborting\n");
        abort();
    }
}

/* evaluate_testdata_mtch: check that testdata read mapped to expected target */
int eval_testdata_mtch(testdata_vars * testdata, char * td_trg_id, int crnt_mut_cd, char * end5p, h_node **p_rdnd)
{    
    if (crnt_mut_cd == NAT) {  //current read contains native sequence
        testdata->nat_mpd++;   //increment native mapped counter
        
        //check that the id of the target used to map the read
        //matches that of the target used to generate the read
        if (!strcmp(td_trg_id, (*p_rdnd)->trg->id)) {
            testdata->nat_exp++;  //id match, increment native expected target counter
        } else {
            testdata->nat_err++;  //mismatch, increment native error counter
            
            //print error describing id mismatch
            printf("TESTDATA: ERROR - native sequence read mapped to incorrect target.\n");
            printf("READ id:  %s\n", td_trg_id);
            printf("TRGT id:  %s\n", (*p_rdnd)->trg->id);
            printf("READ seq: %s\n", end5p);
            printf("TRGT seq: %s\n\n", (*p_rdnd)->trg->sq);
            
            return 0;
        }
        
    } else if (crnt_mut_cd == SUB || //read contains sub/ins/del/mutation
               crnt_mut_cd == INS ||
               crnt_mut_cd == DEL) {
        
        
        //print error describing mapped mutant read
        if (!strcmp(td_trg_id, (*p_rdnd)->trg->id)) {
            printf("TESTDATA: ERROR - mutant (code:%d; sub=0/ins=1/del=2) sequence read mapped to source target.\n", crnt_mut_cd);
            printf("READ id:  %s\n", td_trg_id);
            printf("TRGT id:  %s\n", (*p_rdnd)->trg->id);
            printf("READ seq: %s\n", end5p);
            printf("TRGT seq: %s\n\n", (*p_rdnd)->trg->sq);
            testdata->mut_mpd++; //increment mutant mapped counter
            
            return 0;
            
        } else {
            printf("TESTDATA: NOTE - mutant (code:%d; sub=0/ins=1/del=2) sequence read mapped to non-source target.\n", crnt_mut_cd);
            printf("          This can happen if two variant templates have highly similar sequences\n");
            printf("          The read will not be counted toward the matches for the source target\n");
            printf("READ id:  %s\n", td_trg_id);
            printf("TRGT id:  %s\n", (*p_rdnd)->trg->id);
            printf("READ seq: %s\n", end5p);
            printf("TRGT seq: %s\n\n", (*p_rdnd)->trg->sq);
            //mut_mpd is not incremented because non-source target matches are ignored
            
            return -1; //return -1 can be used to subtract matched read from target matches count
        }
        
        
    } else {
        printf("map_expected_reads: error - unrecognized test data mutation code. aborting...\n");
        abort();
    }
    
    return 0; //this isn't reachable
}

/* print_testdata_analysis: assess testdata analysis
 outcomes and print results to a file and to the screen */
void print_testdata_analysis(TDSPLY_metrics *met, testdata_vars * testdata, target_params * trg_prms, target * trgts)
{
    int i = 0;                //general purpose index
    char out_str[8192] = {0}; //string for storing output line
    
    char pass[7] = "[PASS]";  //output string for PASS outcomes
    char fail[7] = "[FAIL]";  //output string for FAIL outcomes
    char * pf_ptr = NULL;     //pointer to PASS/FAIL strings
    
    //open testdata analysis output file
    FILE * testdata_fp = NULL;
    if ((testdata_fp = fopen("test_data_analysis.txt", "w")) == NULL) {
        printf("print_testdata_analysis: ERROR - could not generate test data analysis file. aborting ...\n");
        abort();
    }
    
    printf("\n");
    sprintf(out_str, "***** TEST DATA ANALYSIS SUMMARY *****\n");
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report number targets assessed
    sprintf(out_str, "> targets - %d targets were assessed\n", trg_prms->t_cnt);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    sprintf(out_str, "> targets - %d non-redundant targets were used for test data analysis\n", trg_prms->nr_cnt);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    sprintf(out_str, "> targets - %d redundant targets were ignored\n", trg_prms->rdndnt);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report number of processed test data reads
    sprintf(out_str, ">   depth - %d total sequences were analyzed\n", met->reads_processed);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    sprintf(out_str, ">   depth - %d native sequences were analyzed\n", testdata->nat_rds);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    sprintf(out_str, ">   depth - %d mutant sequences were analyzed\n", testdata->mut_rds);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    //report mapping efficiency
    pf_ptr = (testdata->nat_mpd == testdata->nat_rds) ? &pass[0] : &fail[0];
    sprintf(out_str, "> mapping - %s %10.6f%% (expected 100%%) of native reads mapped to a target\n",
            pf_ptr, ((float)(testdata->nat_mpd)/(float)(testdata->nat_rds))*100);
    printf2_scrn_n_fl(testdata_fp, out_str);
            
    pf_ptr = (testdata->nat_exp == testdata->nat_rds) ? &pass[0] : &fail[0];
    sprintf(out_str, "> mapping - %s %10.6f%% (expected 100%%) of native reads mapped to the expected target\n",
            pf_ptr, ((float)(testdata->nat_exp)/(float)(testdata->nat_rds))*100);
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    pf_ptr = (testdata->mut_mpd == 0) ? &pass[0] : &fail[0];
    sprintf(out_str, "> mapping - %s %10.6f%% (expected 0.0%%) of mutant reads mapped to a target\n",
            pf_ptr, ((float)(testdata->mut_mpd)/(float)(testdata->mut_rds))*100);
    printf2_scrn_n_fl(testdata_fp, out_str);
              
    //report target coverage
    for (i = 0, pf_ptr = &pass[0]; i < trg_prms->t_cnt; i++) {
        //check whether reads mapped to each target
        //targets with mul=FALSE were the first instance of non-unique targets
        //and were therefore entered into the hash table and used to map reads
        if (!trgts[i].mul && !trgts[i].cnt) {
            sprintf(out_str, "> mapping - [WARNING] zero sequences mapped to target %s\n", trgts[i].id);
            printf2_scrn_n_fl(testdata_fp, out_str);
            pf_ptr = &fail[0];
        }
    }
    
    if (pf_ptr[1] == 'P') {
        sprintf(out_str, "> mapping - %s every target was used when mapping sequences\n", pf_ptr);
    } else {
        sprintf(out_str, "> mapping - %s not every target was used when mapping sequences\n", pf_ptr);
    }
    printf2_scrn_n_fl(testdata_fp, out_str);
    
    
    double test_val = 0; //used to store calculated values
    double expctd_chan[3] = {0.4, 0.4, 0.2}; //expected channel distribution for native reads
    char chan_nm[3][32] = {"bound", "unbound", "unmappable"}; //channel names
    
    for (i = 0; i < CHANNEL_MAX; i++) {
        //calculate fraction of reads in channel i
        test_val = (double)(met->chan_count[i])/(double)(met->matches);
        
        //check if fraction of reads in channel i matches expected value
        pf_ptr = (compare_float(test_val, expctd_chan[i], 0.000001)) ? &pass[0] : &fail[0];
        
        //print output
        sprintf(out_str, "> channel - %s %10.6f%% (expected  %f%%) of mapped sequences contained a(n) %s barcode\n",
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
    
    //close testdata analysis output file
    if ((fclose(testdata_fp)) == EOF) {
        printf("print_3pEnd_testdata_analysis: error - error occurred when closing test data analysis file. Aborting program...\n");
        abort();
    }
}

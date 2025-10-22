//
//  assess_TDSPLY_test_data.c
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
#include "../../seq_utils/mapping_metrics.h"

#include "assess_TDSPLY_test_data.h"

/* parse_testdata_id: parse test data id line and store attributes */
int parse_testdata_id(testdata_vars * testdata, char **td_trg_id, int * crnt_mut_cd, char * id_line)
{
    extern const char refeq[5]; //string that indicates reference id
    extern const char vareq[5]; //string that indicates variant id
    
    int i = 0; //general purpose index
    
    char * ref_loc = NULL;    //pointer to 'REF=' string in test data id line
    char * ref_indx = NULL;   //pointer to reference index in test data id line
    char * var_loc = NULL;    //pointer to 'VAR=' string in test data id line
    char * td_mut_cd = NULL;  //pointer to mutation code in test data id line
    
    //find 'REF=' string
    if ((ref_loc = strstr(id_line, refeq)) == NULL) { //locate refeq identifier in test data read id
        printf("parse_testdata_id: error - unrecognized testdata read id format. aborting\n");
        abort();
    }
    ref_indx = &ref_loc[strlen(refeq)]; //set start of ref index, ref index is tested below to confirm it is all digits
    
    
    //find 'VAR=' string
    if ((var_loc = strstr(id_line, vareq)) == NULL) { //locate vareq identifier in test data read id
        printf("parse_testdata_id: error - unrecognized testdata read id format. aborting\n");
        abort();
    }
    
    //set start of variant id
    if (isdigit(var_loc[strlen(vareq)]) || var_loc[strlen(vareq)] == '-') {
        *td_trg_id = &var_loc[strlen(vareq)];
    } else {
        printf("parse_testdata_id: error - unrecognized testdata read id format. aborting\n");
        abort();
    }
    
    //terminate ref index string
    if (((long long unsigned int)ref_indx < (long long unsigned int)var_loc) && (var_loc[-1] == '_')) {
        var_loc[-1] = '\0';
    } else {
        printf("parse_testdata_id: error - unrecognized testdata read id format. aborting\n");
        abort();
    }
    
    //check that ref index string is composed of digits
    for (i = 0; ref_indx[i]; i++) {
        if (!isdigit(ref_indx[i])) {
            printf("parse_testdata_id: error - unrecognized testdata read id format. aborting\n");
            abort();
        }
    }
    
    
    //set mutation code
    if ((td_mut_cd = strstr(*td_trg_id, "NAT")) != NULL) {        //test data read contains native sequence
        testdata->nat_rds++;                                   //increment native test data read count
        *crnt_mut_cd = NAT;                                    //set current mutation code to native
        
    } else if ((td_mut_cd = strstr(*td_trg_id, "SUB")) != NULL) { //test data read is a substitution mutant
        testdata->mut_rds++;                                   //increment mutant test data read count
        *crnt_mut_cd = SUB;                                    //set current mutation code to substitution
        
    } else if ((td_mut_cd = strstr(*td_trg_id, "INS")) != NULL) { //test data read is an insertion mutant
        testdata->mut_rds++;                                   //increment mutant test data read count
        *crnt_mut_cd = INS;                                    //set current mutation code to insertion
        
    } else if ((td_mut_cd = strstr(*td_trg_id, "DEL")) != NULL) { //test data read is a deletion mutant
        testdata->mut_rds++;                                   //increment mutant test data read count
        *crnt_mut_cd = DEL;                                    //set current mutation code to deletion
        
    } else {  //unrecognized mutation code
        printf("parse_testdata_id: error - could not find mutation code in read id. aborting\n");
        abort();
    }
    
    //set end of target id
    if (td_mut_cd[-1] == '_') { //mutation code should be preceded by a '_'
        td_mut_cd[-1] = '\0';   //terminate target id string
    } else {
        printf("parse_testdata_id: error - unrecognized testdata read id format. aborting\n");
        abort();
    }

    return atoi(ref_indx);
}

/* evaluate_testdata_mtch: check that testdata read mapped to expected target */
int eval_testdata_mtch(testdata_vars * testdata, int td_ref_indx, char * td_trg_id, int crnt_mut_cd, char * end5p, h_node **p_rdnd, target * refs)
{
    int i = 0;       //general purpose index
    int match = 1;   //match flag
    int ret_val = 0; //return value
    
    opt_mx_trg * trg_val = (opt_mx_trg *)(*p_rdnd)->trg->opt;
    target * mapd_ref = trg_val->ref;
    
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
            print_testdata_idNsq_msg(td_trg_id, end5p, (*p_rdnd)->trg->id, (*p_rdnd)->trg->sq);
            
            ret_val = 0;
        }
        
    } else if (crnt_mut_cd == SUB || //read contains sub/ins/del/mutation
               crnt_mut_cd == INS ||
               crnt_mut_cd == DEL) {
        
        
        //print error describing mapped mutant read
        if (!strcmp(td_trg_id, (*p_rdnd)->trg->id)) {
            printf("TESTDATA: ERROR - mutant (code:%d; sub=0/ins=1/del=2) sequence read mapped to source target.\n", crnt_mut_cd);
            
            //duplication of the terminal base by an insertion causes the insertion mutant to be mappable.
            //the code below detects whether this has occurred and, if it has, sets a return value that
            //causes the error to be corrected. if the error occurred for any other reason, no correction
            //is performed.
            
            if (crnt_mut_cd == INS) { //if the current read is an insertion mutant
                
                //test whether the read matches the sequence of the target until
                //the end of one of the two strings is reached.
                for (i = 0, match = 1; end5p[i] && (*p_rdnd)->trg->sq[i] && match; i++) {
                    if (end5p[i] != (*p_rdnd)->trg->sq[i]) { //if a mismatch is found
                        match = 0;                           //set match to false
                    }
                }
                
                if (match) {                                       //if the sequences match
                    if (end5p[i] && !(*p_rdnd)->trg->sq[i]) {      //and the read is longer than the target
                        if (end5p[i] == (*p_rdnd)->trg->sq[i-1]) { //and the 3' nt of the target was duped in the read
                            
                            //print explanation of error and set return val to -1,
                            //which indicates that correction should be perfomed
                            printf("                  This error arose because an insertion at the 3' end of the target duplicated the terminal base.\n");
                            printf("                  The read will not be counted toward the matches for the source target\n");
                            ret_val = -1;
                        }
                    }
                }
            }
            
            print_testdata_idNsq_msg(td_trg_id, end5p, (*p_rdnd)->trg->id, (*p_rdnd)->trg->sq);
            
            if (ret_val != -1) {     //flag to perform correction was not set
                testdata->mut_mpd++; //increment mutant mapped counter
                ret_val = 0;         //set return to zero (no correction to be performed)
            }
            
        } else if ((((long long unsigned int)(&refs[td_ref_indx])) != ((long long unsigned int)mapd_ref))) {
            printf("TESTDATA: ERROR - mutant (code:%d; sub=0/ins=1/del=2) sequence read mapped to non-source target.\n", crnt_mut_cd);
            printf("          This happened because two variant templates have highly similar sequences\n");
            printf("          and a testdata read that contains a mutation mapped to a non-source target\n");
            printf("          The read will not be counted toward the matches for the source target\n");
            print_testdata_idNsq_msg(td_trg_id, end5p, (*p_rdnd)->trg->id, (*p_rdnd)->trg->sq);
            //mut_mpd is not incremented because non-source target matches are ignored
            
            ret_val = -1; //return -1 can be used to subtract matched read from target matches count
            
        } else {
            printf("TESTDATA: ERROR - unrecognized error case. please contact estrobel@buffalo.edu so we can figure out the error.\n");
            print_testdata_idNsq_msg(td_trg_id, end5p, (*p_rdnd)->trg->id, (*p_rdnd)->trg->sq);
        }
        
        
    } else {
        printf("parse_testdata_id: error - unrecognized test data mutation code. aborting...\n");
        abort();
    }
    
    return ret_val; //this isn't reachable
}

/* print_testdata_idNsq_msg: prints id and sequence info for testdata error messges */
void print_testdata_idNsq_msg(char * rd_id, char * rd_sq, char * trg_id, char * trg_sq)
{
    printf("READ id:  %s\n", rd_id);    //print read id
    printf("TRGT id:  %s\n", trg_id);   //print target id
    printf("READ seq: %s\n", rd_sq);    //print read seq
    printf("TRGT seq: %s\n\n", trg_sq); //print target seq
}

/* print_testdata_analysis: assess testdata analysis
 outcomes and print results to a file and to the screen */
void print_testdata_analysis(mapping_metrics *met, testdata_vars * testdata, target_params * trg_prms, target * trgts)
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
    
    for (i = 0; i < TDSPLY_CHANNEL_MAX; i++) {
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

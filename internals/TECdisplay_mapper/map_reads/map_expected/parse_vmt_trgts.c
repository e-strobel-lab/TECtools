//
//  parse_vmt_trgts.c
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

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/revcomp.h"
#include "../../../seq_utils/isDNAbase.h"
#include "../../../seq_utils/isIUPACbase.h"
#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"
#include "../../../seq_utils/basemap.h"

#include "./set_trgt.h"
#include "../../../cotrans_preprocessor/MUX_trgt_gen/set_barcoded_compact_target.h"

#include "parse_vmt_trgts.h"

/* parse_header_lines: parse targets file header lines for expected variant count
 and wild type sequence information */
void parse_header_lines(FILE * ifp, target_params * trg_prms, TDSPLY_fasta * wt)
{
    int i = 0; //general purpose index
    int j = 0; //general pursose index
    
    char line1[MAX_LINE] = {0}; //array for storing line 1 of targets file
    char line2[MAX_LINE] = {0}; //array for storing line 2 of targets file
    
    char * xpctd_trgs = NULL;   //pointer to the start of expected targets count string
    char * wt_nm = NULL;        //pointer to the start of wild type sequence name
    char * wt_sq = NULL;        //pointer to the start of wild type sequence
    
    const char xpctd_code[10] = "variants:"; //variant count code
    const char wt_code[4] = "WT:";           //wild type sequence code
    
    //get expected variant count from targets file line 1
    get_line(line1, ifp); //get first line of targets file, which contains expected number of targets
    if (!memcmp(line1, xpctd_code, strlen(xpctd_code)) && line1[strlen(xpctd_code)]) { //check line1 formatting
        xpctd_trgs = &line1[strlen(xpctd_code)];          //set pointer to start of expected target count
        for (i = 0; xpctd_trgs[i] && i < MAX_LINE; i++) { //check that all xpctd_trgs chars are digits
            if (!isdigit(xpctd_trgs[i])) {
                printf("parse_header_lines: error - expected targets count must be a string of digits. aborting...");
                abort();
            }
        }
        trg_prms->xpctd = atoi(xpctd_trgs); //set the expected number of targets
        
    } else {
        printf("parse_header_lines: error - unexpected format for variant count line. aborting...");
        abort();
    }
    
    //get wt sequence information from targets file line 2
    get_line(line2, ifp); //get second line of targets file, which contains name and wt sequence
    if (!memcmp(line2, wt_code, strlen(wt_code)) && line2[strlen(wt_code)]) {
        wt_nm = &line2[strlen(wt_code)];
        for (i = 0; wt_nm[i] != '\t' && wt_nm[i] && i < MAX_LINE; i++) {;} //iterate to tab delimiter
        
        if (!line2[i] || i == MAX_LINE) { //check loop exit condition for unexpected line format
            printf("parse_header_lines: error - unexpected wt seq line format (missing tab delimiter). aborting...\n");
            abort();
        } else {
            wt_nm[i] = '\0';     //replace tab with null to split the line
            wt_sq = &wt_nm[i+1]; //set pointer to target sequence (starts at index i+1)
        }
        
        //allocate memory for wt name
        if ((wt->nm = malloc((strlen(wt_nm)+1) * sizeof(*(wt->nm)))) == NULL) {
            printf("parse_header_lines: error - memory allocation for wt name failed. aborting...\n");
            abort();
        }
        
        //allocate memory for wt seq
        if ((wt->sq = malloc((strlen(wt_sq)+1) * sizeof(*(wt->sq)))) == NULL) {
            printf("parse_header_lines: error - memory allocation for wt sequence failed. aborting...\n");
            abort();
        }
        
        strcpy(wt->nm, wt_nm); //store wt name
        strcpy(wt->sq, wt_sq); //store wt seq
        
    } else {
        printf("parse_header_lines: error - unexpected format for wt sequence line. aborting...\n");
        abort();
    }
    
    return;
}

/* parse_vmt_trgts: parse targets file to obtain target ids, sequences, attributes,
  and min/max transcript lengths */
void parse_vmt_trgts(FILE * ifp, int trgt_ftype, target * refs, opt_ref * ref_val, void * trgts, void * trg_val, target_params * trg_prms, TDSPLY_fasta * wt, int mode)
{
    //TODO: double check that reversion from revcomp is all correct
    extern int debug;				           //flag to run debug mode
    
    int i = 0;                                 //general purpose index
    int j = 0;                                 //general purpose index
    
    char line[MAX_LINE];			           //array for input line
    int got_line = 0;                          //flag to check success of get_line function
    
    char *trgt_id = NULL;			           //pointer to target name in target id
    char *trgt_sq = NULL;			           //pointer to target sequence in target id
    char unprcsd_ref_sq[MAXLEN+1] = {0};       //array to store unprocessed reference sequence
    
    target * crnt_ref = NULL;                  //pointer to current reference target,
    opt_ref * crnt_ref_val = NULL;             //pointer to current reference target values
    opt_mx_trg * crnt_trg_val = NULL;          //pointer to current target values
    
    const char ref_code[5] = "REF:";           //reference target code
    const char tpr_code[5] = "TPR:";           //targets per reference code
    const char vbs_code[5] = "VBS:";           //variable bases code
    const char cnst_code[7] = "const:";        //constants code
    
    char * tmp_ptr = NULL;                     //temporary pointer to reference name component
    char * ref_name_ptr = NULL;                //pointer to reference name in reference target id
    char * xpctd_tpr_ptr = NULL;               //pointer to expected targets per reference string in reference target id
    char * vbases_ptr = NULL;                  //pointer to variable bases string in reference target id
    char * cnstnts_ptr = NULL;                 //pointer to string of constant indels in reference target id
    int xpctd_tpr = 0;                         //expected targets per reference
    int crrnt_tpr = 0;                         //number of targets that have been processed for current reference target
      
    printf("***********************  target parsing  ************************\n");
    
    while (get_line(line, ifp)) {
        
        //if t_cnt equals xpctd and got_line was successful,
        //there are more targets than expected based on the
        //contents of targets file line 1 and not enough
        //memory was allocated.
        if (trg_prms->t_cnt > trg_prms->xpctd) {
            printf("parse_vmt_trgts: error - more targets than expected in targets file. aborting...\n");
            abort();
        }
                
        //split input line into target name and sequence
        for (i = 0; line[i] != '\t' && line[i] && i < MAX_LINE; i++) { ;} //iterate to first tab
        if (!line[i] || i == MAX_LINE) { //check loop exit condition for unexpected line format
            printf("parse_vmt_trgts: error - unexpected end target line format (missing tab delimiter). aborting...\n");
            abort();
        } else {
            line[i] = '\0';	      //replace tab with null to split line
            trgt_id = &line[0];   //set pointer to target id (starts at index 0)
            trgt_sq = &line[i+1]; //set pointer to target sequence (starts at index i+1)
        }
                        
        if (!memcmp(trgt_id, ref_code, strlen(ref_code)) && trgt_id[strlen(ref_code)]) { //check if ref target
        
            //check that correct number of targets were processed for previous reference sequence
            check_tpr_match(trg_prms->r_cnt, trg_prms->t_per_r[trg_prms->r_cnt-1], xpctd_tpr);
            
            if (trg_prms->r_cnt >= MAXREF) { //check that maximum reference target numer has not been reached
                printf("parse_vmt_trgts: error - more reference targets than expected maximum (%d) in targets file. aborting...\n", MAXREF);
                abort();
            }
            
            //parse reference target id. '|' characters separate reference target attributes
            ref_name_ptr = &trgt_id[strlen(ref_code)]; //set pointer to start of reference target name
            tmp_ptr = NULL;                            //set tmp_ptr to NULL
            xpctd_tpr_ptr = NULL;                      //set expected tpr pointer to NULL
            vbases_ptr = NULL;                         //set vbases pointer to NULL
            cnstnts_ptr = NULL;                        //set constants pointer to NULL
            
            for (i = 0; trgt_id[i]; i++) {
                if (trgt_id[i] == '|') {          //found separator
                    trgt_id[i] = '\0';            //split target id string at separator
                    if (trgt_id[i+1]) {           //separator is followed by an attribute
                        tmp_ptr = &trgt_id[++i];  //set pointer to the start of the attribute
                        
                        //determine attribute type.
                        //transcripts per reference attributes are preceeded by "TPR:"
                        //constant indels attributes are preceeded by "const:"
                        
                        if (!memcmp(tmp_ptr, tpr_code, strlen(tpr_code)) && tmp_ptr[strlen(tpr_code)]) {
                            xpctd_tpr_ptr = &tmp_ptr[strlen(tpr_code)]; //set pointer to TPR string
                        
                        } else if (!memcmp(tmp_ptr, vbs_code, strlen(vbs_code)) && tmp_ptr[strlen(vbs_code)]) {
                            vbases_ptr = &tmp_ptr[strlen(vbs_code)];    //set pointer to variable bases string
                            
                        } else if (!memcmp(tmp_ptr, cnst_code, strlen(cnst_code)) && tmp_ptr[strlen(cnst_code)]) {
                            cnstnts_ptr = &tmp_ptr[strlen(cnst_code)];  //set pointer to constant indels string
                            
                        } else { //unrecognized target id format
                            printf("parse_vmt_trgts: error - unrecognized target id format for reference target %s\n. aborting...\n", ref_name_ptr);
                            abort();
                        }
                        
                    } else { //target is missing attribute
                        printf("parse_vmt_trgts: error - separator character in target id for reference target %s is not followed by an attribute. aborting...\n", ref_name_ptr);
                        abort();
                    }
                }
            }
                                    
            //check that tpr is composed of digits and
            //convert the expected tpr string to an int
            for (i = 0; xpctd_tpr_ptr[i]; i++) {
                if (!isdigit(xpctd_tpr_ptr[i])) {
                    printf("parse_vmt_trgts: error - transcripts per reference value for reference target %s contains a non-digit character. aborting...\n", ref_name_ptr);
                    abort();
                }
            }
            xpctd_tpr = atoi(xpctd_tpr_ptr); //set the number of expected targets for the current reference
            crrnt_tpr = 0;                   //initialize the current targets per reference counter
            strcpy(unprcsd_ref_sq, trgt_sq); //store unprocessed reference sequence for target seq validation
            
            validate_ref_seq(ref_name_ptr, trgt_sq, wt); //check that reference sequence is valid
            
            //set target structure values for the current reference target
            set_ref(&refs[trg_prms->r_cnt], &ref_val[trg_prms->r_cnt], ref_name_ptr, trgt_sq, xpctd_tpr, vbases_ptr, cnstnts_ptr);
            crnt_ref = &refs[trg_prms->r_cnt]; //set crnt_ref pointer to point to current reference target struct
            if (debug) {print_reference_debug(&refs[trg_prms->r_cnt], trg_prms);} //print debug messages
            
            trg_prms->r_cnt++; //increment the reference target counter
                        
        } else if (trg_prms->r_cnt) { //check that a reference target has been processed
            
            //confirm target is a match to reference seq
            //this confirms that the target matches its reference at all non-variable
            //positions and that target sequence length matches reference sequencel length
            if (!validate_trgt_seq(trgt_id, trgt_sq, unprcsd_ref_sq)) {
                printf("parse_vmt_trgts: error - target %s is associated with reference %s but is not a complete match at non-variable positions. aborting...\n", trgt_id, crnt_ref->id);
                abort();
            }
            
            //set target structure values for the current target
            if (mode == TDSPLY_TRGS) {
                set_trgt(&(((target *)trgts)[trg_prms->t_cnt]),
                         &(((opt_mx_trg *)trg_val)[trg_prms->t_cnt]),
                         crnt_ref, trgt_id, trgt_sq);
                if (debug) {print_target_debug(&(((target *)trgts)[trg_prms->t_cnt]), trg_prms);} //print debug messages
            } else if (mode == TPROBE_TRGS) {
                set_barcoded_compact_target(&(((compact_target *)trgts)[trg_prms->t_cnt]),
                                            &(((opt_BC *)trg_val)[trg_prms->t_cnt]),
                                            crnt_ref, trgt_id, trgt_sq, trg_prms, trgt_ftype);
            }
            
            trg_prms->t_cnt++;                       //increment parsed target counter
            trg_prms->t_per_r[trg_prms->r_cnt-1]++;  //increment targets per reference counter
            
        } else { //file must begin with a reference target. throw error and abort.
            printf("parse_vmt_trgts: error - first target must be a reference target. aborting...\n");
            abort();
        }
    }
    
    //check that correct number of targets were processed for previous reference sequence
    check_tpr_match(trg_prms->r_cnt, trg_prms->t_per_r[trg_prms->r_cnt-1], xpctd_tpr);
    
    //if t_cnt is less than xpctd, the targets file
    //did not contain the correct number of targets
    if (trg_prms->t_cnt != trg_prms->xpctd) { //t_cnt should be equal to xpctd
        printf("parse_vmt_trgts: error - incorrect number (%d) of targets in targets file. expected (%d). aborting...\n", trg_prms->t_cnt, trg_prms->xpctd);
        abort();
    }

    printf("the targets file contained %d sequences associated with %d reference targets\n", trg_prms->t_cnt, trg_prms->r_cnt);
    
    return;
}

/* check_tpr_match: check that expected tpr matches actual tpr */
void check_tpr_match(int cnt, int actual, int xpctd)
{
    if (actual != xpctd) {
        printf("check_tpr_match: error - actual number of targets for reference %d (%d) is not equal to the expected number indicated in targets file (%d). aborting...", cnt, actual, xpctd);
        abort();
    } else if (cnt) {
        printf("reference target %d: %d/%d\t(%5.2f%%) of expected targets were processed\n", cnt, actual, xpctd, ((float)(actual)/(float)(xpctd))*100);
    }
    
    return;
}

/* validate_ref_seq: check that reference comprises valid characters and is not too long */
void validate_ref_seq(char * ref_nm, char * ref_sq, TDSPLY_fasta * wt)
{
    int i = 0; //general purpose index
    
    //check reference target length
    if (strlen(ref_sq) > MAXLEN) {
        printf("validate_ref_seq: error - reference target %s is %lu nts long. expected <=%d nts. aborting...\n", ref_nm, strlen(ref_sq), MAXLEN);
    }
    
    //check that reference target length matches wt sequence length
    if (strlen(ref_sq) != strlen(wt->sq)) {
        printf("validate_ref_seq: error - length of reference target %s (%lu nts) is not equal to the length of the aligned wt sequence (%lu nts). aborting...\n", ref_nm, strlen(ref_sq), strlen(wt->sq));
    }
    
    //check that reference target seq comprises valid characters
    for (i = 0; ref_sq[i]; i++) {
        if (!isIUPACbase(ref_sq[i]) && ref_sq[i] != '.' && ref_sq[i] != '-') {
            printf("validate_ref_seq: error - reference %s contains invalid character %c. aborting...\n", ref_nm, ref_sq[i]);
            abort();
        }
    }
    
    return;
}


/* validate_trgt_seq: check that target comprises valid characters and
 matches its associated reference target at all non-variable positions */
int validate_trgt_seq(char * trgt_nm, char * trgt_sq, char * ref_sq)
{
    int i = 0; //general purpose index
 
    //check that target length matches reference target length
    if (strlen(trgt_sq) != strlen(ref_sq)) {
        return 0; //return fail
    }
    
    for (i = 0; trgt_sq[i] && ref_sq[i] && i < MAXLEN; i++) {
        
        //check that target seq comprises valid characters
        if (!isDNAbase(trgt_sq[i]) && trgt_sq[i] != '.' && trgt_sq[i] != '-') {
            printf("validate_trgt_seq: error - target %s contains invalid character %c. aborting...\n", trgt_nm, trgt_sq[i]);
            abort();
        }
        
        //check that target matches reference target at non-variable positions
        if (isDNAbase(ref_sq[i]) || ref_sq[i] == '.' || ref_sq[i] == '-') {
            if (toupper(trgt_sq[i]) != toupper(ref_sq[i])) {
                return 0; //return fail
            }
        }
    }
    
    return 1; //return success
}

/* process_trgt_seq: convert target seq to upper case and remove non-base characters. */
int process_trgt_seq(char * ipt, char * processed_seq)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
        
    for (i = 0, j = 0; ipt[i] && i < MAXLEN; i++) {
        if (isIUPACbase(ipt[i])) {                //if the input char is a IUPAC base,
            processed_seq[j++] = toupper(ipt[i]); //copy it to the filtered sequence as uppercase
        }
    }
    
    processed_seq[i] = '\0'; //terminate the filtered sequence string
        
    return 1; //return success
}

/* set_ref: set reference target values in targets structure */
void set_ref(target * refs, opt_ref * ref_val, char * ref_id, char * ref_sq, int tpr, char * vbases, char * cnstnts)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    opt_ref * p_ref_val = NULL;	//pointer for dereferencing opt pointer in targets structure to opt_ref
    
    char processed_seq[MAXLEN+1] = {0}; //array to store seq w/o non-base chars
    
    if (!process_trgt_seq(ref_sq, &processed_seq[0])) { //filter non-base chars from input seq
        printf("set_ref: error - non-base filtering failed. aborting...\n");
        abort();
    }
    
    //allocate memory for target structure id and sequence values
    if ((refs->id = malloc((strlen(ref_id)+1) * sizeof(*(refs->id)))) == NULL) {
        printf("set_ref: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    if ((refs->sq = malloc((strlen(processed_seq)+1) * sizeof(*(refs->sq)))) == NULL) {
        printf("set_ref: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    //set reference id and seq values target structure
    strcpy(refs->id, ref_id);	      //copy reference id to target entry
    strcpy(refs->sq, processed_seq);  //copy processed reference sequence to target entry
    refs->opt = ref_val;		      //set refs_opt pointer to address of current opt_ref struct
    p_ref_val = (opt_ref *)refs->opt; //dereference refs->opt to simplify code below
    
    //set transcripts per reference value
    p_ref_val->tpr = tpr;
    
    //set unprocessed input sequence
    if ((p_ref_val->ipt_sq = malloc((strlen(ref_sq)+1) * sizeof(*(p_ref_val->ipt_sq)))) == NULL) {
        printf("set_ref: error - memory allocation failed. aborting...\n");
        abort();
    }
    strcpy(p_ref_val->ipt_sq, ref_sq);
    
    //set variable bases string
    if ((p_ref_val->vbases = malloc((strlen(vbases)+1) * sizeof(*(p_ref_val->vbases)))) == NULL) {
        printf("set_ref: error - memory allocation failed. aborting...\n");
        abort();
    }
    strcpy(p_ref_val->vbases, vbases);
    
    //if present, set constant indels string
    if (cnstnts != NULL) {
        if ((p_ref_val->cnstnts = malloc((strlen(cnstnts)+1) * sizeof(*(p_ref_val->cnstnts)))) == NULL) {
            printf("set_ref: error - memory allocation failed. aborting...\n");
            abort();
        }
        strcpy(p_ref_val->cnstnts, cnstnts);
    }
    
    /* copy indices of variable bases in reference
     sequence to the current ref_val struct */
    for (i = 0; refs->sq[i] && i < MAX_LINE; i++) {
        if (!isDNAbase(refs->sq[i])) {       //non-native DNA base indicates variable position
            if (!isIUPACbase(refs->sq[i])) { //variable base must be an IUPAC base
                printf("set_ref: error - unrecognized base %c in reference target sequence. aborting\n", refs->sq[i]);
                abort();
            }
    
            //set the position of variable bases in the reference sequence
            if (p_ref_val->vb_cnt < SEQ2BIN_MAX_KEY) {      //check that number of variable bases does not exceed max
                p_ref_val->vb_pos[p_ref_val->vb_cnt++] = i; //set variable base position in reference target
            } else { //if number of variable bases exceeds max, throw error and abort
                printf("set_ref: number of variable bases (%d) exceeds allowed maximum (%d). aborting...\n", p_ref_val->vb_cnt, SEQ2BIN_MAX_KEY);
                abort();
            }
        }
    }
    
    return;
}

/* print_reference_debug: print debug messages for reference processing outcome */
void print_reference_debug(target * trgt, target_params * trg_prms)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    opt_ref * crnt_ref_val = (opt_ref *)trgt->opt; //dereference reference target values
    
    printf("ref target %d\n", trg_prms->r_cnt);    //print reference target number
    printf("%s\n", trgt->id);                      //print reference target id
    printf("TPR:\t%d\n", crnt_ref_val->tpr);       //print transcripts per reference
    
    if (crnt_ref_val->cnstnts != NULL) {           //if present, print constants string
        printf("const:\t%s\n", crnt_ref_val->cnstnts);
    }
    
    for (i = 0; i < crnt_ref_val->vb_cnt; i++) {   //print variable base positions
        printf("%3d ", crnt_ref_val->vb_pos[i]);
    }
    
    printf("\n%s\n", trgt->sq); //print reference target sequence
    
    //print variable bases aligned to target sequence
    for (i = 0, j = 0; trgt->sq[i]; i++) {   //for every base of the reference sequence
        if (i == crnt_ref_val->vb_pos[j]) {  //if the position is a variable base position
            printf("%c", trgt->sq[crnt_ref_val->vb_pos[j++]]); //print the identity of the base, increment j
        } else {
            printf("."); //otherwise print a '.'
        }
    }
    printf("\n\n");
        
    return;
}

/* print_target_debug: print debug messages for target processing outcome */
void print_target_debug(target * trgt, target_params * trg_prms)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    opt_mx_trg * crnt_trg_val = NULL; //pointer to target values
    opt_ref * crnt_ref_val = NULL;    //pointer to reference target values
    
    crnt_trg_val = (opt_mx_trg *)trgt->opt;                //set pointer to target values
    crnt_ref_val = (opt_ref *)crnt_trg_val->ref->opt;      //set pointer to reference target values
    
    printf("target %d\t", trg_prms->t_cnt);                //print target number
    printf("%s\t%s\n%s\n", trgt->id, trgt->key, trgt->sq); //print target id and sequence
    
    //print key aligned to target sequence
    for (i = 0, j = 0; trgt->sq[i]; i++) {   //for every base of the target sequence
        if (i == crnt_ref_val->vb_pos[j]) {  //if the position is a variable base position
            printf("%c", trgt->key[j++]);    //print the identity of the base
        } else {
            printf(".");                     //otherwise print a '.'
        }
    }
    printf("\n\n");

    return;
}

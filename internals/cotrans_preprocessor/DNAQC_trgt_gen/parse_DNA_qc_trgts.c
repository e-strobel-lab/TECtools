//
//  parse_DNA_qc_trgts.c
//  
//
//  Created by Eric Strobel on 9/8/26.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"
#include "../../variant_maker/bc_ind.h"

#include "../MUX_trgt_gen/set_barcoded_compact_target.h"

#include "./link_trgts.h"

#include "parse_DNA_qc_trgts.h"

/* parse_DNA_qc_trgts: parse DNA qc targets fasta file */
void parse_DNA_qc_trgts(FILE * ifp, int trgt_ftype, compact_target * ctrg, opt_BC * trg_val, association * assoc, target_params * trg_prms)
{
    extern char bc1_ind[3]; //barcode 1 indicator
    extern char bc2_ind[4]; //barcode 2 indicator
    
    uint64_t i = 0; //general purpose index
    uint64_t a = 0; //association index
    
    compact_target * p_trgt1 = NULL; //pointer to compact target 1 of a pair
    compact_target * p_trgt2 = NULL; //pointer to compact target 2 of a pair
        
    char trgt_id[MAX_LINE+1] = {0};        //array for current target identifier
    char trgt_sq[MAX_LINE+1] = {0};        //current target sequence
    char prsd_bc_sq[MAX_LINE+1] = {0};     //target sequence string to be parsed for barcodes
    char prsd_insert_sq[MAX_LINE+1] = {0}; //target sequence string to be parsed for the insert
    
    char * var_id = NULL; //pointer to variant id
    char * bc_id1 = NULL; //pointer to barcode 1 id
    char * bc_id2 = NULL; //pointer to barcode 2 id
    
    char * bc1 = NULL; //pointer to barcode 1
    char * bc2 = NULL; //pointer to barcode 2
    
    char * insert = NULL; //pointer to target insert
        
    int L1 = 0; //variable for checking get fasta line 1 success
    int L2 = 0; //variable for checking get fasta line 2 success
    
    int brcd_len; //barcode length
    
    while ((L1 = get_line(trgt_id, ifp))) { //get fasta target id line
        
        //if t_cnt > xpctd and got_line was successful,
        //there are more targets than expected and not
        //enough memory was allocated.
        if (trg_prms->t_cnt > trg_prms->xpctd) {
            printf("parse_DNA_qc_trgts: error - more targets than expected in targets file. aborting...\n");
            abort();
        }
        
        L2 = get_line(trgt_sq, ifp);     //get fasta sequence line
        strcpy(prsd_bc_sq, trgt_sq);     //make copy of sequence line for barcode parsing
        strcpy(prsd_insert_sq, trgt_sq); //make copy of sequence line for insert parsing
        
        if (L1 == 0 || L2 == 0) { //if getting line 1 or line 2 failed, throw error and abort
            printf("parse_DNA_qc_trgts: error - reached end of barcodes file before reading expected number of barcodes. aborting...\n");
            abort();
        }
        
        parse_DNA_qc_trgt_id(&var_id, &bc_id1, &bc_id2, trgt_id); //parse target id
        brcd_len = parse_DNA_qc_barcodes(&bc1, &bc2, prsd_bc_sq); //parse barcodes
        parse_DNA_qc_insert(&insert, prsd_insert_sq, brcd_len);   //parse insert
        
        if (!trg_prms->BClen) {         //if the barcode length was not yet set
            trg_prms->BClen = brcd_len; //set brcd_len
        }
        
        /* set_DNA_qc_compact_target: set compact_target values for DNA qc target */
        void set_DNA_qc_compact_target(compact_target * ctrg, opt_BC * BC_val, char * insert, char * var_id, char * bc_id, char * bc_sq, int bc_num, target_params * trg_prms);
        
        //set barcode 1 compact_target
        p_trgt1 = &ctrg[trg_prms->t_cnt]; //set pointer to barcode 1 compact_target for use below
        set_DNA_qc_compact_target(p_trgt1, &trg_val[trg_prms->t_cnt], insert, var_id, bc_id1, bc1, 1, trg_prms);
        trg_prms->t_cnt++; //increment target count
        
        //set barcode 2 compact_target
        p_trgt2 = &ctrg[trg_prms->t_cnt]; //set pointer to barcode 2 compact_target for use below
        set_DNA_qc_compact_target(p_trgt2, &trg_val[trg_prms->t_cnt], insert, var_id, bc_id2, bc2, 2, trg_prms);
        trg_prms->t_cnt++; //increment target count
        
        link_trgts(p_trgt1, p_trgt2, LINK_NTVS); //link targets 1 and 2 as a pair
        
        //set pointers to association structures
        p_trgt1->utl = &assoc[a++];
        p_trgt2->utl = &assoc[a++];
                
        L1 = L2 = 0;       //zero L1 and L2
    }
    
    return;
}

/* parse_DNA_qc_trgt_id: parse target id of DNA qc target */
void parse_DNA_qc_trgt_id(char ** var_id, char ** bc_id1, char ** bc_id2, char * trgt_id)
{
    extern char bc1_ind[3]; //barcode 1 indicator
    extern char bc2_ind[4]; //barcode 2 indicator
    
    int i = 0; //general purpose index
    
    //set start of variant id to the character after the leading '>'
    if (trgt_id[0] == '>') {
        (*var_id) = &trgt_id[1];
    } else {
        printf("parse_DNA_qc_trgt_id: error - expected target sequences to be in fasta format. aborting...\n");
        abort();
    }
    
    //iterate to the first '_', which precedes the first barcode id,
    //and terminate the variant id string
    for (i = 0; trgt_id[i] && trgt_id[i] != '_'; i++) {;}
    if (trgt_id[i] == '_') {
        trgt_id[i++] = '\0';
    } else {
        printf("parse_DNA_qc_trgt_id: error - unexpected target id format. aborting...\n");
        abort();
    }
        
    //if the barcode 1 indicator follows the underscore identified above,
    //set pointer to the numerical value of barcode 1 id.
    if (!memcmp(&trgt_id[i], bc1_ind, strlen(bc1_ind))) {
        (*bc_id1) = &trgt_id[i+strlen(bc1_ind)];
    } else {
        printf("parse_DNA_qc_trgt_id: error - unexpected target id format. aborting...\n");
        abort();
    }
    
    //iterate to next '_', which should be preceded by digits, and set terminating null
    for (i += strlen(bc1_ind); trgt_id[i] && isdigit(trgt_id[i]); i++) {;}
    if (trgt_id[i] == '_') {
        trgt_id[i++] = '\0';
    } else {
        printf("parse_DNA_qc_trgt_id: error - unexpected target id format. aborting...\n");
        abort();
    }
    
    //if the barcode 2 indicator follows the underscore identified above,
    //set pointer to the numerical value of barcode 2 id.
    if (!memcmp(&trgt_id[i], bc2_ind, strlen(bc2_ind))) {
        (*bc_id2) = &trgt_id[i+strlen(bc2_ind)];
    } else {
        printf("parse_DNA_qc_trgt_id: error - unexpected target id format. aborting...\n");
        abort();
    }
    
    //iterate to the end of the string, which should be preceded by digits.
    //throw error if the end of the string is not reached
    for (i += strlen(bc2_ind); trgt_id[i] && isdigit(trgt_id[i]); i++) {;}
    if (trgt_id[i]) {
        printf("parse_DNA_qc_trgts: error - unexpected format for target id. aborting...\n");
        abort();
    }
    
    return;
}

/* parse_DNA_qc_barcodes: parse barcode sequences of DNA qc target */
int parse_DNA_qc_barcodes(char ** bc1, char ** bc2, char * prsd_bc_sq)
{
    extern char vra5[22];
    extern char RLA29synch_3p11[34];
    
    char * p_vra5  = NULL; //pointer to vra5 sequence
    char * p_synch = NULL; //pointer to synchronization sequence
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    //identify vra5 and point bc1 to the character after the vra5 sequence
    if ((p_vra5 = strstr(prsd_bc_sq, vra5)) != NULL) {
        *bc1 = &p_vra5[strlen(vra5)];
    } else {
        printf("parse_DNA_qc_barcodes: error - did not detect vra5 sequence in target sequence. aborting...\n");
        abort();
    }
    
    //identify the synchronization site and point bc2 to the character after the site
    if ((p_synch = strstr(prsd_bc_sq, RLA29synch_3p11)) != NULL) {
        *bc2 = &p_synch[strlen(RLA29synch_3p11)];
    } else {
        printf("parse_DNA_qc_barcodes: error - did not detect synchronization site in target sequence. aborting...\n");
        abort();
    }
    
    //both barcodes are expected to be 16 nucleotides long (MAX_BARCODE_LEN)
    //iterate until a lower case character is reached to set minimum barcode length
    for (i = 0; isupper((*bc1)[i]) && (*bc1)[i]; i++) {;}
    for (j = 0; isupper((*bc2)[j]) && (*bc2)[j]; j++) {;}
    
    //if the barcodes are the expected length, set terminating null and return the barcode length
    if (i == MAX_BARCODE_LEN && j == MAX_BARCODE_LEN) {
        (*bc1)[MAX_BARCODE_LEN] = '\0';
        (*bc2)[MAX_BARCODE_LEN] = '\0';
        return i;
        
    } else { //otherwise, throw error and abort
        printf("parse_DNA_qc_barcodes: error - expected barcodes to be %d nucleotides long. aborting...\n", MAX_BARCODE_LEN);
        abort();
    }
}

/* parse_DNA_qc_insert: parse insert of DNA qc target */
void parse_DNA_qc_insert(char **insert, char * prsd_insert_sq, int brcd_len)
{
    extern char vra5[22];
    extern char RLA29synch_3p11[34];
    
    char * p_vra5  = NULL; //pointer to vra5 sequence
    char * p_synch = NULL; //pointer to synchronization sequence
    
    //set pointer to the start of the insert
    if ((p_vra5 = strstr(prsd_insert_sq, vra5)) != NULL) {
        *insert = &p_vra5[strlen(vra5) + brcd_len];
    } else {
        printf("parse_DNA_qc_insert: error - did not detect vra5 sequence in target sequence. aborting...\n");
        abort();
    }
    
    //set character after the insert to terminating null
    if ((p_synch = strstr(prsd_insert_sq, RLA29synch_3p11)) != NULL) {
        p_synch[strlen(RLA29synch_3p11)] = '\0';
    } else {
        printf("parse_DNA_qc_insert: error - did not detect synchronization site in target sequence. aborting...\n");
        abort();
    }
    
    return;
}

/* set_DNA_qc_compact_target: set compact_target values for DNA qc target */
void set_DNA_qc_compact_target(compact_target * ctrg, opt_BC * BC_val, char * insert, char * var_id, char * bc_id, char * bc_sq, int bc_num, target_params * trg_prms)
{
    extern char bc1_ind[3]; //barcode 1 indicator
    extern char bc2_ind[4]; //barcode 2 indicator
    
    int i = 0; //general purpose index
    
    //check that barcode sequence is composed of uppercase DNA bases
    for (i = 0; bc_sq[i]; i++) {
        if (!isDNAbase(bc_sq[i]) || !isupper(bc_sq[i])) {
            printf("set_DNA_qc_compact_target: barcode '%s' contains a formatting error. aborting...\n", insert);
            abort();
        }
    }
    
    char * p_bc_ind = NULL; //barcode indicator pointer

    if (bc_num == 1) {        //if setting barcode 1 target
        p_bc_ind = bc1_ind;   //set pointer to barcode 1 indicator
    } else if (bc_num == 2) { //if setting barcode 2 target
        p_bc_ind = bc2_ind;   //set pointer to barcode 2 indicator
    } else {
        printf("set_DNA_qc_compact_target: error - barcode indicator must be '1' or '2'. aborting...\n");
        abort();
    }
    
    uint64_t nid = 0;           //numerical barcode id
    char fid[MAX_LINE+1] = {0}; //full id
    int ret = 0;                //snprintf return value
    
    //generate full barcode id
    ret = snprintf(fid, MAX_LINE, "%s_%s%s", var_id, p_bc_ind, bc_id);
    if (ret >= MAX_LINE || ret < 0) {
        printf("set_DNA_qc_compact_target: error - failed to generate target full id. aborting...\n");
        abort();
    }
    
    set_BC_val(ctrg, BC_val, insert, NULL, NATIVE_BRCD); //set pointer to BC_val structure
    nid = (uint64_t)(strtoull(bc_id, NULL, 10));         //convert nid string to numerical id
    ctrg->bid = nid << MUTCODE_BITS;     //set barcode id (bits 0-15 are used for mutCode)
    seq2bin_long(bc_sq, &ctrg->bsq, 1);  //store barcode sequence as a binary sequence
    BC_val->typ = NAT;                   //set barcode type to native
    BC_val->num = bc_num;                //set barcode number
    
    //allocate memory for storing character-encoded barcode
    if ((ctrg->csq = malloc((strlen(bc_sq)+1) * sizeof(*ctrg->csq))) == NULL) {
        printf("set_barcoded_compact_target: error - memory allocation for character-encoded target sequence failed. aborting...\n");
        abort();
    }
    strcpy(ctrg->csq, bc_sq); //store character-encoded barcode

    //allocate memory for storing full barcode id
    if ((ctrg->cid = malloc((strlen(fid)+1) * sizeof(*ctrg->cid))) == NULL) {
        printf("set_barcoded_compact_target: error - memory allocation for full barcode id failed. aborting...\n");
        abort();
    }
    strcpy(ctrg->cid, fid); //store character-encoded barcode id
    
    return;    
}

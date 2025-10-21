//
//  set_barcoded_compact_target.c
//  
//
//  Created by Eric Strobel on 10/1/25.
//

#include <stdio.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"

#include "../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"

#include "set_barcoded_compact_target.h"


/* set_barcoded_compact_target: set target values in target struct */
void set_barcoded_compact_target(compact_target * ctrg, opt_BC * BC_val, target * crnt_ref, char * trgt_id, char * trgt_sq, target_params * trg_prms, int trgt_ftype, int data_type)
{
    extern char pra1_sc1[41];
    extern char c3sc1[34];
    extern char vra3[27];
    extern char RLA29synch_3p11[34]; //v2
    
    char * p_nid = NULL;  //pointer to numerical barcode id
    char * p_fid = NULL;  //pointer to full barcode id
    char * p_trgt = NULL; //pointer to target RNA sequence
    char * p_brcd = NULL; //pointer to barcode sequence
    
    uint64_t nid = 0; //numerical barcode id
    
    char * ldr2use = NULL;  //leader sequence to search for when parsing target
    char * trlr2use = NULL; //trailer sequence to search for when parsing target
    char * lnkr2use = NULL; //linker sequence to search for when parsing target
    
    if (data_type == TDSPLY) { //when processing TECdisplay data
        ldr2use = c3sc1;       //leader sequence is C3-SC1
        trlr2use = vra3;       //trailer sequence is VRA3
        lnkr2use = NULL;       //no linker is used //TODO: add ability to search for user-specified linker?
        
    } else if (data_type == TPROBE_MUX) { //when processing TECprobe-MUX data
        ldr2use = pra1_sc1;               //leader sequence is PRA1-SC1 adapter
        trlr2use = vra3;                  //trailer is VRA3
        lnkr2use = RLA29synch_3p11;       //linker is RLA29synch_3p11
    }
    
    parse_barcode_id(&p_nid, &p_fid, trgt_id, trgt_ftype); //parse barcode id to isolate full and numerical ids
    parse_MUX_target_seq(&p_trgt, &p_brcd, &trg_prms->BClen, trgt_sq, trgt_ftype, ldr2use, trlr2use, lnkr2use); //parse seq to isolate target and BC seqs
    
    //printf("%s\t%s\t%s\t%s\n", p_nid, p_fid, p_brcd, p_trgt);
    
    set_BC_val(ctrg, BC_val, p_trgt, NULL, NATIVE_BRCD); //set pointer to BC_val structure
    nid = (uint64_t)(strtoull(p_nid, NULL, 10));         //convert nid string to numerical id
    ctrg->bid = nid << MUTCODE_BITS;     //set barcode id (bits 0-15 are used for mutCode)
    seq2bin_long(p_brcd, &ctrg->bsq, 1); //store barcode sequence as a binary sequence
    BC_val->typ = NAT;                   //set barcode type to native
    
    //allocate memory for storing character-encoded barcode
    if ((ctrg->csq = malloc((strlen(p_brcd)+1) * sizeof(*ctrg->csq))) == NULL) {
        printf("set_barcoded_compact_target: error - memory allocation for character-encoded target sequence failed. aborting...\n");
        abort();
    }
    strcpy(ctrg->csq, p_brcd); //store character-encoded barcode

    //allocate memory for storing full barcode id (only performed for input barcodes)
    if ((ctrg->cid = malloc((strlen(p_fid)+1) * sizeof(*ctrg->cid))) == NULL) {
        printf("set_barcoded_compact_target: error - memory allocation for full barcode id failed. aborting...\n");
        abort();
    }
    strcpy(ctrg->cid, p_fid); //store character-encoded barcode id
    
    return;
}

/* parse_barcode_id: parse barcode id to set full id and numerical id pointers */
int parse_barcode_id(char ** p_nid, char ** p_fid, char * crnt_bcid, int trgt_ftype) {
    
    uint64_t i = 0; //general purpose index
    
    char * p_last_uscore = NULL; //pointer to last detected underscore
    int only_digits = 0;         //flag that barcode is composed entirely of digits
    
    if (trgt_ftype == VMT_FILE) {   //processing VMT file
        if (crnt_bcid[0] != '>') {  //check that target id does not start with '>'
            *p_fid = &crnt_bcid[0]; //set pointer to full barcode id
        } else { //if fasta entry starts with '>', throw error and abort
            printf("parse_barcode_id: error - unexpected format for vmt file. was a fasta file provided instead? aborting...\n");
            abort();
        }
        
    } else if (trgt_ftype == FASTA_FILE) { //processing fasta file
        if (crnt_bcid[0] == '>') {         //check that target id begins with '>'
            *p_fid = &crnt_bcid[1];        //set pointer to char after '>'
        } else {
            printf("parse_barcode_id: error - unexpected format for fasta file. was a vmt file provided instead? aborting...\n");
            abort();
        }
        
    } else {
        printf("parse_barcode_id: error - unrecognized targets file type. aborting...\n");
        abort();
    }
    
    //get numerical barcode id by searching for last underscore and
    //setting pointer to the digit that follows the last underscore
    for (i = 0, p_last_uscore = NULL, only_digits = 1; crnt_bcid[i]; i++) {
        
        if (crnt_bcid[i] == '_') {            //if underscore is found
            p_last_uscore = &crnt_bcid[i];    //set pointer to most recent underscore
            only_digits = 1;                  //reset only digits flag to true
            
        } else if (!isdigit(crnt_bcid[i])) {  //if a non-digit character was found
            only_digits = 0;                  //set only digits flag to false
        }
    }
    
    if (i > 0 && p_last_uscore == NULL && only_digits) { //id is strictly numerical
        *p_nid = *p_fid;     //set numerical id pointer to id start of full id
        return NUMERICAL_ID; //return NUMERICAL_ID code
        
    } else if (i > 0 && p_last_uscore != NULL && only_digits) { //id contains more than just a number
        *p_nid = &p_last_uscore[1]; //set numerical id pointer to one character after the last underscore
        return COMPLEX_ID;          //return COMPLEX_ID code
        
    } else { //unexpected id format, throw error and abort
        printf("parse_barcode_id: unexpected id format. aborting...\n");
        abort();
    }
}

/* parse_MUX_target_seq: parse TECprobe-MUX target sequence to identify barcode and target RNA sequences */
void parse_MUX_target_seq(char ** p_trgt, char ** p_brcd, int * brcd_len, char * crnt_seq, int trgt_ftype, char * ldr, char * trlr, char * lnkr)
{
    char * p_ldr = NULL;  //pointer to leader
    char * p_trlr = NULL; //pointer to trailer
    char * p_lnkr = NULL; //pointer to linker
    
    char tmp_brcd[MAX_BARCODE_LEN+1] = {0}; //array to temporarily store barcode during transposition
    
    if (trgt_ftype == VMT_FILE) { //vmt file sequences contain only the target, linker, and barcode
        *p_trgt = &crnt_seq[0];   //set target start to first array member
        
    } else if (trgt_ftype == FASTA_FILE) { //fasta file sequences also contain leader and trailer sequences
        
        //set pointer to leader sequence
        if ((p_ldr = strstr(crnt_seq, ldr)) == NULL) {
            printf("parse_target_seq: leader not detected in target sequence in fasta file mode. aborting...\n");
            abort();
        }
        
        //check that leader sequence is at the start of the target sequence
        if ((uint64_t)(&p_ldr[0]) == (uint64_t)(&crnt_seq[0])) {  //if pointer is at the start of the target seq
            *p_trgt = &crnt_seq[strlen(ldr)];                     //set pointer to first nt after leader
        } else {
            printf("parse_target_seq: leader not detected at head of target sequence in fasta file mode. aborting...\n");
            abort();
        }
        
        //set pointer to trailer sequence
        if ((p_trlr = strstr(*p_trgt, trlr)) == NULL) {
            printf("parse_target_seq: trailer not detected in target sequence in fasta file mode. aborting...\n");
            abort();
        }
        p_trlr[0] = '\0'; //terminate barcode string by zeroing first index of trailer string, which follows barcode
        
    } else {
        printf("parse_target_seq: error - unrecognized target file type. aborting...");
        abort();
    }
    
    
    //set pointer to linker
    if (lnkr != NULL) {
        if ((p_lnkr = strstr(*p_trgt, lnkr)) == NULL) {
            printf("parse_target_seq: linker not detected. aborting...\n");
            abort();
        }
        p_lnkr[0] = '\0'; //terminate target string by zeroing first index of linker string
        *p_brcd = &p_lnkr[strlen(lnkr)];   //set barcode pointer to nt that follows linker
        (*p_brcd)[strlen(*p_brcd)] = '\0'; //guarantee null char after barcode
        
    } else {
        if (strlen(*p_trgt) > MAX_BARCODE_LEN) {                       //check that target is longer than BC
            *p_brcd = &((*p_trgt)[strlen(*p_trgt) - MAX_BARCODE_LEN]); //set pointer to start of barcode
            strcpy(tmp_brcd, *p_brcd);                                 //store temp copy of barcode
            (*p_brcd)[0] = '\0';                                       //terminate target string
            *p_brcd = &((*p_brcd)[1]);                                 //point BC to char after trg str null char
            strcpy(*p_brcd, tmp_brcd);                                 //copy barcode back to seq string
            (*p_brcd)[strlen(*p_brcd)] = '\0';                         //guarantee null char after barcode
            
        } else {
            printf("parse_target_seq: unexpected short target. aborting...\n");
            abort();
        }
    }
    
    //set or check barcode length
    if (!(*brcd_len)) {                     //if barcode length was not previously set,
        *brcd_len = strlen(*p_brcd);        //set barcode length
        if (*brcd_len != MAX_BARCODE_LEN) { //check that barcode length matches MAX_BARCODE_LEN
            printf("parse_target_seq: unexpected barcode length (%d nt, expected %d nt). aborting...\n", *brcd_len, MAX_BARCODE_LEN);
            abort();
        }
        
        //if linker was provided, check that barcode length matches
        //the distance between the linker and the trailer sequences
        if (trgt_ftype == FASTA_FILE && lnkr != NULL) {
            if (*brcd_len != (uint64_t)(&p_trlr[0]) - (uint64_t)(&(*p_brcd)[0])) {
                printf("parse_target_seq: error - unexpected distance betwen end of linker and start of trailer. aborting...\n");
                abort();
            }
        }
        
    } else if (strlen(*p_brcd) != *brcd_len) { //otherwise, confirm current BC len matches first BC len
        printf("parse_target_seq: discordant barcode lengths detected. aborting...\n");
        abort();
    }
    
    return;
}

/* set_BC_val: set optional values for barcode target */
void set_BC_val(compact_target * ctrg, opt_BC * BC_val, char * tsq, compact_target * ntv, int mode)
{
    opt_BC * p_opt_BC = NULL; //pointer for handling barcode target optional values
    
    ctrg->opt = BC_val;             //point opt to corresponding barcode values
    p_opt_BC = (opt_BC *)ctrg->opt; //set pointer for handling BC_val below
        
    char processed_seq[MAX_LINE+1] = {0}; //array for storing processed target sequence
    
    int len = 0; //target sequence length
    
    if (mode == NATIVE_BRCD) { //if setting values for native barcode
        p_opt_BC->ref = ctrg;  //point reference to self
        len = strlen(tsq);     //determine target sequence length
        
        //convert target seq to uppercase and remove spacing characters
        if (len > MAX_LINE) {
            printf("set_BC_val: error - target sequence length exceeds array bounds. aborting...\n");
            abort();
        } else {
            process_trgt_seq(tsq, processed_seq);
        }
        
        //allocate memory for target sequence and store target sequence in barcode target optional values
        if ((p_opt_BC->tsq = malloc((len+1) * sizeof(*p_opt_BC->tsq))) == NULL) {
            printf("set_BC_val: error - memory allocation for target sequence failed. aborting...\n");
            abort();
        }
        strcpy(p_opt_BC->tsq, processed_seq);
        
    } else if (mode == MUTANT_BRCD) { //if setting values for mutant barcode
        p_opt_BC->ref = ntv;          //point ref to reference native barcode from which the mutant was derived
    }
    
    return;
}

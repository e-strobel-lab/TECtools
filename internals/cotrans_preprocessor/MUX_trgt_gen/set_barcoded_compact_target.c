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
void set_barcoded_compact_target(compact_target * ctrg, opt_BC * BC_val, target * crnt_ref, char * trgt_id, char * trgt_sq, target_params * trg_prms, int trgt_ftype)
{
    char * p_nid = NULL;  //pointer to numerical barcode id
    char * p_fid = NULL;  //pointer to full barcode id
    char * p_trgt = NULL; //pointer to target RNA sequence
    char * p_brcd = NULL; //pointer to barcode sequence
    
    uint64_t nid = 0; //numerical barcode id
    
    parse_barcode_id(&p_nid, &p_fid, trgt_id, trgt_ftype); //parse barcode id to isolate full and numerical ids
    parse_MUX_target_seq(&p_trgt, &p_brcd, &trg_prms->BClen, trgt_sq, trgt_ftype); //parse seq to isolate target and BC seqs
    
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
void parse_MUX_target_seq(char ** p_trgt, char ** p_brcd, int * brcd_len, char * crnt_seq, int trgt_ftype)
{
    extern char pra1_sc1[41];        //PRA1_SC1 adapter
    extern char vra3[27];            //VRA3_adapter
    extern char RLA29synch_3p11[34]; //RLA29synch_3p11 linker
    
    char * p_pra1_adptr = NULL; //pointer to PRA1_SC1 adapter
    char * p_vra3 = NULL;       //pointer to VRA3 adapter
    char * p_lnkr = NULL;       //pointer to linker
    
    if (trgt_ftype == VMT_FILE) { //vmt file sequences contain only the target, linker, and barcode
        *p_trgt = &crnt_seq[0];   //set target start to first array member
        
    } else if (trgt_ftype == FASTA_FILE) { //fasta file sequences also contain PRA1-SC1 and VRA3 sequences
        
        //set pointer to PRA1_SC1 adapter sequence
        if ((p_pra1_adptr = strstr(crnt_seq, pra1_sc1)) == NULL) {
            printf("parse_target_seq: PRA1_SC1 adapter not detected in target sequence in fasta file mode. aborting...\n");
            abort();
        }
        
        //check that PRA1_SC1 sequence is at the start of the oligo sequence
        if ((uint64_t)(&p_pra1_adptr[0]) == (uint64_t)(&crnt_seq[0])) {  //if pointer is at the start of the oligo seq
            *p_trgt = &crnt_seq[strlen(pra1_sc1)];                       //set pointer to first nt after PRA1_SC1
        } else {
            printf("parse_target_seq: PRA1_SC1 sequence not detected at head of target sequence in fasta file mode. aborting...\n");
            abort();
        }
        
        //set pointer to VRA3 adapter sequence
        if ((p_vra3 = strstr(crnt_seq, vra3)) == NULL) {
            printf("parse_target_seq: VRA3 adapter not detected in target sequence in fasta file mode. aborting...\n");
            abort();
        }
        p_vra3[0] = '\0'; //terminate barcode string by zeroing first index of VRA3 adapter string, which follows barcode
        
    } else {
        printf("parse_target_seq: error - unrecognized target file type. aborting...");
        abort();
    }
    
    //set pointer to RLA29synch_3p11 linker
    if ((p_lnkr = strstr(crnt_seq, RLA29synch_3p11)) == NULL) {
        printf("parse_target_seq: RLA29synch_3p11 linker not detected. aborting...\n");
        abort();
    }
    p_lnkr[0] = '\0'; //terminate target string by zeroing first index of linker string, which follows target string
    *p_brcd = &p_lnkr[strlen(RLA29synch_3p11)]; //set barcode pointer to nt that follows RLA29synch_3p11 linker
    
    if (!(*brcd_len)) {
        *brcd_len = strlen(*p_brcd); //get barcode length
        if (*brcd_len != MAX_BARCODE_LEN) {
            printf("parse_target_seq: unexpected barcode length (%d nt, expected %d nt). aborting...\n", *brcd_len, MAX_BARCODE_LEN);
            abort();
        }
        
        if (trgt_ftype == FASTA_FILE) {
            //check that barcode length matches distance between the RLA29synch_3p11 linker and VRA3 adapter
            if (*brcd_len != (uint64_t)(&p_vra3[0]) - (uint64_t)(&(*p_brcd)[0])) {
                printf("parse_target_seq: error - unexpected distance betwen end of RLA29synch_3p11 linker and start of VRA3 adapter. aborting...\n");
                abort();
            }
        }
        
    } else if (strlen(*p_brcd) != *brcd_len) {
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

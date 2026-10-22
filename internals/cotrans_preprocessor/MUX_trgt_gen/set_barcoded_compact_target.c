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
#include "../../seq_utils/C3SCbc_tree.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"

#include "../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"

#include "set_barcoded_compact_target.h"


/* set_barcoded_compact_target: set target values in target struct */
void set_barcoded_compact_target(compact_target * ctrg, opt_BC * BC_val, target * crnt_ref, char * trgt_id, char * trgt_sq, target_params * trg_prms, int trgt_ftype, int lib_struct, int data_type)
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
    
    if (data_type == TDSPLY_STD || data_type == TDSPLY_BRCD) { //when processing TECdisplay data
        ldr2use  = c3sc1;           //leader sequence is C3-SC1
        trlr2use = vra3;            //trailer sequence is VRA3
        lnkr2use = NULL;            //no linker is used //TODO: add ability to search for user-specified linker?
        
    } else if (data_type == TPROBE_MUX && lib_struct == SGL_BC) { //when processing single BC TECprobe-MUX data
        ldr2use  = pra1_sc1;        //leader sequence is PRA1-SC1 adapter
        trlr2use = vra3;            //trailer is VRA3
        lnkr2use = RLA29synch_3p11; //linker is RLA29synch_3p11
        
    } else if (data_type == TPROBE_MUX && lib_struct == DBL_BC) {
        ldr2use  = NULL;            //leader is processed by parse_C3SC_barcode
        trlr2use = vra3;            //trailer is VRA3
        lnkr2use = RLA29synch_3p11; //linker is RLA29synch_3p11
        
    } else {
        printf("set_barcoded_compact_target: error - unrecongnized data type/library structure. aborting...\n");
        abort();
    }
    
    static C3SC_barcode_node scbt_root = {EOF, NULL, NULL, NULL, NULL}; //root for sc barcode tree
    C3SC_barcode_node * scbc_node = NULL;                               //pointer for storing sc barcode node return
    
    if (lib_struct == DBL_BC) {                    //if processing doubly-barcoded targets
        if (scbt_root.key == EOF) {                //if the tree has not yet been generated
            make_C3SCbc_tree(&scbt_root, 0, 0);    //generate the tree
        }
        
        //parse the sc barcode and find match in sc barcode tree
        if ((scbc_node = parse_C3SC_barcode(&trgt_sq, &scbt_root)) == NULL) {
            printf("set_barcoded_compact_target: error - failed to find match for SC barcode during target parsing. aborting...\n");
            abort();
        }
    }
    
    parse_barcode_id(&p_nid, &p_fid, trgt_id, trgt_ftype); //parse barcode id to isolate full and numerical ids
    parse_MUX_target_seq(&p_trgt, &p_brcd, &trg_prms->BClen, trgt_sq, trgt_ftype, ldr2use, trlr2use, lnkr2use); //parse seq to isolate target and BC seqs
    
    //printf("%s\t%s\t%s\t%s\n", p_nid, p_fid, p_brcd, p_trgt);
    
    set_BC_val(ctrg, BC_val, p_trgt, NULL, lib_struct, scbc_node, NATIVE_BRCD); //set pointer to BC_val structure
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

/* parse_C3SC_barcode: parse C3SC barcode and find match in sc barcode tree. sc barcode
parsing also points trgt_sq to the start of the target sequence */
C3SC_barcode_node * parse_C3SC_barcode(char ** trgt_sq, C3SC_barcode_node * scbt_root)
{
    extern char c3[20];
    extern char ttcg[5];
    extern char c3sc1[34];
    
    char * trgt_cpy = NULL;   //copy of target sequence
    char * p_c3 = NULL;       //pointer to C3 sequence
    char * p_ustrm_bc = NULL; //pointer to upstream barcode segment
    char * p_ttcg = NULL;     //pointer to ttcg loop sequence
    char * p_dstrm_bc = NULL; //pointer to downstream barcode sequence
    char * p_trgt = NULL;     //pointer to target sequence
    
    C3SC_barcode_node * scbc_node = NULL; //pointer for storing outcome of sc barcode tree search
    
    //make copy of target sequence for parsing
    if ((trgt_cpy = malloc((strlen(*trgt_sq)+1) * sizeof(*trgt_cpy))) == NULL) {
        printf("parse_C3SC_barcode: failed to allocate memory for target copy. aborting...\n");
        abort();
    }
    strcpy(trgt_cpy, *trgt_sq);
    
    //trim c3 sequence from input sequence
    p_c3 = strstr(*trgt_sq, c3);                                    //find C3 sequence in target
    if (p_c3 == NULL || (uint64_t)(p_c3) != (uint64_t)(*trgt_sq)) { //verify C3 sequence formatting
        printf("parse_C3SC_barcode: error - failed to detect C3 sequence\n. ");
        abort();
    }
    *trgt_sq = &p_c3[strlen(c3)];
    
    //parse sc barcode from target sequence copy
    p_c3 = strstr(trgt_cpy, c3);                                                 //find C3 sequence in target
    p_trgt     = &(trgt_cpy[strlen(c3sc1)]);                                     //set ptr to start of target seq
    p_ustrm_bc = &(trgt_cpy[strlen(c3)]);                                        //set ptr to upstream bc stem
    p_ttcg     = &(trgt_cpy[strlen(c3) + MAX_C3SCBC_TREE_LAYER]);                //set ptr to ttcg loop
    p_dstrm_bc = &(trgt_cpy[strlen(c3) + MAX_C3SCBC_TREE_LAYER + strlen(ttcg)]); //set ptr to downstream bc stem
    
    //verify upstream bounds of target sequence, which should be preceded by two 'a' nucleotides
    if (p_trgt[-1] == 'a' && p_trgt[-2] == 'a') { //if upstream bounds are verified
        p_trgt[-2] = '\0';                        //terminate downstream barcode string
        
    } else { //otherwise, throw error and abort
        printf("parse_C3SC_barcode: error - unexpected leader structure. 'aa' motif not detected before target start. aborting...\n");
        abort();
    }
        
    //verify ttcg tetraloop
    if (!memcmp(p_ttcg, ttcg, strlen(ttcg))) { //if ttcg tetraloop is verified
        p_ttcg[0] = '\0';                      //terminate upstream barcode string
        
    } else { //otherwise, throw error and abort
        printf("parse_C3SC_barcode: error - unexpected leader structure. 'ttcg' not detected between stems of structure cassette. aborting...\n");
        abort();
    }
    
    //perform search for barcode match
    scbc_node = srch_C3SCbc_tree(p_ustrm_bc, 0, &p_ustrm_bc[0], scbt_root);
        
    //verify match node and revcomp match node sequences are as expected
    if (strcmp(scbc_node->sq, p_ustrm_bc) || strcmp(scbc_node->rc->sq, p_dstrm_bc)) {
        printf("parse_C3SC_barcode: error - reported barcode match is not a match. this should never happen. aborting...\n");
        abort();
    }
    
    //free allocated memory
    free(trgt_cpy);
    trgt_cpy = NULL;
    
    return scbc_node; //return node
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
    
    static int ldr_detected = -1; //flag that a leader was or was not detected in the first target
        
    if (trgt_ftype == VMT_FILE) { //vmt file sequences contain only the target, linker, and barcode
        *p_trgt = &crnt_seq[0];   //set target start to first array member
        
    } else if (trgt_ftype == FASTA_FILE) { //fasta file sequences may contain leader and trailer sequences
        
        //the code below is designed to allow flexibility in the provided targets. each assay has an expected
        //leader sequence, and if detected, this leader will be trimmed from the target. however the presence
        //of the leader is not required and non-standard leaders can be used (but will not be trimmed).
        
        //check for presence of leader sequence
        if (ldr != NULL) {                 //if a leader might be present
            p_ldr = strstr(crnt_seq, ldr); //search for leader within the current target sequence
            
            if (ldr_detected == -1) { //if the leader status of the first target has not yet been set
                
                //if a leader was detected in the target
                if (p_ldr != NULL && (uint64_t)(&p_ldr[0]) == (uint64_t)(&crnt_seq[0])) {
                    ldr_detected = 1; //set ldr_detected to true
                } else {              //otherwise,
                    ldr_detected = 0; //set ldr_detected to false
                }
                
            } else if ((p_ldr == NULL && ldr_detected != 0) ||
                       (p_ldr != NULL && ldr_detected == 1)) {
                //if leader status is inconsistent across targets, throw error and abort
                printf("parse_MUX_target_seq: error - inconsistent leader status in targets. aborting...\n");
                abort();
            }
        }
        
        //set pointer to start of target
        if (ldr != NULL && p_ldr != NULL) { //if a leader was expected and identified...
            
            //check that leader sequence is at the start of the target sequence
            if ((uint64_t)(&p_ldr[0]) == (uint64_t)(&crnt_seq[0])) {  //if pointer is at the start of the target seq
                *p_trgt = &crnt_seq[strlen(ldr)];                     //set pointer to first nt after leader
            } else {
                printf("parse_target_seq: leader not detected at head of target sequence in fasta file mode. aborting...\n");
                abort();
            }
            
        } else { //if no leader was identified or expected...
            *p_trgt = &crnt_seq[0]; //set p_trgt to start of crnt_seq
        }
        
        //set pointer to trailer sequence
        if ((p_trlr = strstr(*p_trgt, trlr)) != NULL) {
            p_trlr[0] = '\0'; //terminate barcode string by zeroing first index of trailer string, which follows barcode
        } else {
            printf("parse_MUX_target_seq: trailer sequence was not detected in target sequence. aborting...\n");
            abort();
        }
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
            //^NOTE: all input sequence configurations contain a trailer sequence, which guarantees that there
            //is enough space to shift the barcode down one index in the array
            
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
void set_BC_val(compact_target * ctrg, opt_BC * BC_val, char * tsq, compact_target * ntv, int lib_struct, C3SC_barcode_node * scbc_node, int mode)
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
        
        if (lib_struct == DBL_BC) {
            p_opt_BC->scbc = scbc_node;
        }
        
    } else if (mode == MUTANT_BRCD) { //if setting values for mutant barcode
        p_opt_BC->ref = ntv;          //point ref to reference native barcode from which the mutant was derived
    }
    
    return;
}

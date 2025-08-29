//
//  mk_MUX_trgts.c
//  
//
//  Created by Eric Strobel on 6/25/25.
//

#include <stdio.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"

#include "mk_MUX_trgts.h"

//global variables
uint64_t trgt_ID_mask = 0xFFFFFFFFFFFF0000; //barcode id mask that isolates the reference barcode id
uint64_t mutcode_mask = 0x000000000000FFFF; //barcode id mask that isolates mutation code

/* mk_MUX_trgts: manages TECprobe-MUX target generation */
int mk_MUX_trgts(compact_target * ctrg, opt_BC * BC_val, FILE * fp_MUXtrgs, int brcd_cnt, int clcd_ctrg_cnt)
{
    extern char RLA29synch_3p11[34]; //linker sequence for TECprobe-MUX experiments
    
    int brcd_len = 0;  //barcode length
    int ctrg_cnt = 0;  //number of barcode targets
    
    ctrg_cnt = store_barcode_targets(ctrg, BC_val, &brcd_len, fp_MUXtrgs, brcd_cnt); //store input barcodes as compact targets
    printf("barcode count = %d\ncalc'd trg count = %d\n", brcd_cnt, clcd_ctrg_cnt);
    
    uint64_t i = 0;        //general purpose index
    uint64_t mutCode = 0;  //variable for storing mutation code. each barcode id is an unsigned 64 bit integer.
                           //bits 0-15 are used to store a mutation code. bits >= 16 are used to store the id
                           //of the input barcode from which the mutated barcode is derived
    
    int tot_SUB_trgts = 0; //number of substitution targets generated
    int tot_INS_trgts = 0; //number of insertion targets generated
    int tot_DEL_trgts = 0; //number of deletion targets generated
    
    for (i = 0; i < brcd_cnt; i++) { //for every input barcode
        
        mutCode = 1; //set mutCode to 1 (input barcodes have mutCode 0)
        
        //generate sub and indel targets for the current barcode
        tot_SUB_trgts += mk_SUB_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], brcd_len, &mutCode);
        tot_INS_trgts += mk_INS_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], brcd_len, &mutCode);
        tot_DEL_trgts += mk_DEL_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], brcd_len, &mutCode, toupper(RLA29synch_3p11[strlen(RLA29synch_3p11)-1]));
    }
    
    //print the number of targets that were generated
    printf("total targets = %d\n", brcd_cnt + tot_SUB_trgts + tot_INS_trgts + tot_DEL_trgts);
    printf("%d compact targets\n", ctrg_cnt);
    
    return ctrg_cnt; //return the number of compact targets that were generated
}

/* store_barcode_targets: stores input barcodes as compact targets */
int store_barcode_targets(compact_target * ctrg, opt_BC * BC_val, int * brcd_len, FILE * fp_MUXtrgs, int brcd_cnt)
{
    uint64_t i = 0; //general purpose index
    
    opt_BC * p_opt_BC = NULL; //pointer for handling barcode target optional values
    
    int ctrg_cnt = 0;                  //number of compact targets
    char crnt_bcid[MAX_LINE+1] = {0};  //array for current barcode identifier
    char crnt_oligo[MAX_LINE+1] = {0}; //current oligo sequence
    char line[MAX_LINE+1] = {0};       //array to store lines
    
    char * p_nid = NULL;  //pointer to numerical barcode id
    char * p_fid = NULL;  //pointer to full barcode id
    char * p_trgt = NULL; //pointer to target RNA sequence
    char * p_brcd = NULL; //pointer to barcode sequence
    
    uint64_t nid = 0;     //numerical barcode id
    
    int L1 = 0; //variable for checking get fasta line 1 success
    int L2 = 0; //variable for checking get fasta line 2 success
    
    for (i = 0; i < brcd_cnt; i++) {  //for every barcode in the input barcode file
        
        L1 = get_line(crnt_bcid, fp_MUXtrgs);  //get fasta line 1
        L2 = get_line(crnt_oligo, fp_MUXtrgs); //get fasta line 2
        
        if (L1 == 0 || L2 == 0) { //if getting line 1 or line 2 failed, throw error and abort
            printf("store_barcode_targets: error - reached end of barcodes file before reading expected number of barcodes. aborting...\n");
            abort();
        }
        
        parse_barcode_id(&p_nid, &p_fid, crnt_bcid);             //parse barcode id to isolate full and numerical ids
        parse_oligo_seq(&p_trgt, &p_brcd, brcd_len, crnt_oligo); //parse oligo seq to isolate target and barcode seqs
        
        if (*brcd_len != MAX_BARCODE_LEN) { //if barcode length does not match the max barcode length, throw error and abort
            printf("store_barcode_targets: error - barcode length does not match the maximum barcode length (%d nucleotides), which is currently the only permitted barcode length. aborting...\n", MAX_BARCODE_LEN);
        }
        
        set_BC_val(&ctrg[ctrg_cnt], &BC_val[ctrg_cnt], p_trgt, NULL, NATIVE_BRCD); //set pointer to BC_val structure
        nid = (uint64_t)(strtoull(p_nid, NULL, 10));  //convert nid string to numerical id
        ctrg[ctrg_cnt].bid = nid << MUTCODE_BITS;     //set barcode id (bits 0-15 are used for mutCode)
        seq2bin_long(p_brcd, &ctrg[ctrg_cnt].bsq, 1); //store barcode sequence as a binary sequence
        BC_val[ctrg_cnt].typ = NAT;                   //set barcode type to native
        
        //allocate memory for storing character-encoded barcode
        if ((ctrg[ctrg_cnt].csq = malloc((strlen(p_brcd)+1) * sizeof(*ctrg[ctrg_cnt].csq))) == NULL) {
            printf("store_barcode_targets: error - memory allocation for character-encoded target sequence failed. aborting...\n");
            abort();
        }
        strcpy(ctrg[ctrg_cnt].csq, p_brcd); //store character-encoded barcode
    
        //allocate memory for storing full barcode id (only performed for input barcodes)
        if ((ctrg[ctrg_cnt].cid = malloc((strlen(p_fid)+1) * sizeof(*ctrg[ctrg_cnt].cid))) == NULL) {
            printf("store_barcode_targets: error - memory allocation for full barcode id failed. aborting...\n");
            abort();
        }
        strcpy(ctrg[ctrg_cnt].cid, p_fid); //store character-encoded barcode id
        
        ctrg_cnt++; //increment number of compact targets
    }
    
    //check the end of barcode file was reached, if not, throw error
    if (get_line(line, fp_MUXtrgs)) {
        printf("store_barcode_targets: more barcodes than expected in barcodes file. aborting...\n");
        abort();
    }
    
    return ctrg_cnt; //return number of compact targets that were generated from the input barcodes
}

/* set_BC_val: set optional values for barcode target */
void set_BC_val(compact_target * ctrg, opt_BC * BC_val, char * tsq, compact_target * ntv, int mode)
{
    opt_BC * p_opt_BC = NULL; //pointer for handling barcode target optional values
    
    ctrg->opt = BC_val;             //point opt to corresponding barcode values
    p_opt_BC = (opt_BC *)ctrg->opt; //set pointer for handling BC_val below
        
    int len = 0; //target sequence length
    
    if (mode == NATIVE_BRCD) { //if setting values for native barcode
        p_opt_BC->ref = ctrg;  //point reference to self
        len = strlen(tsq);     //determine target sequence length
        
        //allocate memory for target sequence and store target sequence in barcode target optional values
        if ((p_opt_BC->tsq = malloc((len+1) * sizeof(*p_opt_BC->tsq))) == NULL) {
            printf("set_BC_val: error - memory allocation for target sequence failed. aborting...\n");
            abort();
        }
        strcpy(p_opt_BC->tsq, tsq);
        
    } else if (mode == MUTANT_BRCD) { //if setting values for mutant barcode
        p_opt_BC->ref = ntv;          //point ref to reference native barcode from which the mutant was derived
    }
    
    return;
}

/* parse_barcode_id: parse barcode id to set full id and numerical id pointers */
int parse_barcode_id(char ** p_nid, char ** p_fid, char * crnt_bcid) {
    
    uint64_t i = 0; //general purpose index
    
    char * p_last_uscore = NULL; //pointer to last detected underscore
    int only_digits = 0;         //flag that barcode is composed entirely of digits
    
    if (crnt_bcid[0] == '>') {  //check that fasta entry starts with '>'
        *p_fid = &crnt_bcid[1]; //set pointer to full barcode id
    } else { //if fasta entry does not start with '>', throw error and abort
        printf("store_barcode_targets: error - unexpected format for fasta file. aborting...\n");
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
        *p_nid = &crnt_bcid[1];     //set numerical id pointer to id start
        return NUMERICAL_ID;
        
    } else if (i > 0 && p_last_uscore != NULL && only_digits) { //id contains more than just a number
        *p_nid = &p_last_uscore[1]; //set numerical id pointer to one character after the last underscore
        return COMPLEX_ID;
        
    } else { //unexpected id format, throw error and abort
        printf("store_barcode_targets: unexpected id format. aborting...\n");
        abort();
    }
}

/* parse_oligo_seq: parse oligonucleotide sequence to identify barcode and target RNA sequences */
void parse_oligo_seq(char ** p_trgt, char ** p_brcd, int * brcd_len, char * crnt_oligo)
{
    extern char pra1_sc1[41];        //PRA1_SC1 adapter
    extern char vra3[27];            //VRA3_adapter
    extern char RLA29synch_3p11[34]; //RLA29synch_3p11 linker
    
    char * p_pra1_adptr = NULL; //pointer to PRA1_SC1 adapter
    char * p_vra3 = NULL;       //pointer to VRA3 adapter
    char * p_lnkr = NULL;       //pointer to linker
    
    //set pointer to PRA1_SC1 adapter sequence
    if ((p_pra1_adptr = strstr(crnt_oligo, pra1_sc1)) == NULL) {
        printf("parse_oligo_seq: PRA1_SC1 adapter not detected. aborting...\n");
        abort();
    }
    
    //check that PRA1_SC1 sequence is at the start of the oligo sequence
    if ((uint64_t)(&p_pra1_adptr[0]) == (uint64_t)(&crnt_oligo[0])) {  //if pointer is at the start of the oligo seq
        *p_trgt = &crnt_oligo[strlen(pra1_sc1)];                       //set target pointer to first nt after PRA1_SC1
    } else {
        printf("parse_oligo_seq: PRA1_SC1 sequence not detected at head of oligo sequence. aborting...\n");
        abort();
    }
    
    //set pointer to VRA3 adapter sequence
    if ((p_vra3 = strstr(crnt_oligo, vra3)) == NULL) {
        printf("parse_oligo_seq: VRA3 adapter not detected. aborting...\n");
        abort();
    }
    
    //set pointer to RLA29synch_3p11 linker
    if ((p_lnkr = strstr(crnt_oligo, RLA29synch_3p11)) == NULL) {
        printf("parse_oligo_seq: RLA29synch_3p11 linker not detected. aborting...\n");
        abort();
    }
    *p_brcd = &p_lnkr[strlen(RLA29synch_3p11)]; //set barcode pointer to nt that follows RLA29synch_3p11 linker
    
    p_lnkr[0] = '\0'; //terminate target string by zeroing first index of linker string, which follows target string
    p_vra3[0] = '\0'; //terminate barcode string by zeroing first index of VRA3 adapter string, which follows barcode
    
    *brcd_len = strlen(*p_brcd); //get barcode length
    
    //check that barcode length corresponds to the entire sequence between the RLA29synch_3p11 linker and VRA3 adapter
    if (*brcd_len != (uint64_t)(&p_vra3[0]) - (uint64_t)(&(*p_brcd)[0])) {
        printf("store_barcode_targets: error - unexpected distance betwen end of RLA29synch_3p11 linker and start of VRA3 adapter. aborting...\n");
        abort();
    }
    
    return;
}

/* mk_SUB_trgts: generates single substitution targets for input barcode */
int mk_SUB_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode)
{
    int sub_trgts_made = 0; //number of substitution targets that were made
    
    int pos2var = 0; //position to vary
    int i = 0;       //general purpose index
    
    char var[3][MAX_LINE+1] = {{0}};    //array for constructing substitution targets
    
    if (*mutCode != MIN_SUB_CODE) {
        printf("mk_SUB_trgts: error - incorrect minimum substitution code (%llu). expected %d. aborting...\n", (long long unsigned int)(*mutCode), MIN_SUB_CODE);
        abort();
    }
    
    for (pos2var = 0; pos2var < brcd_len && src_ctrg->csq[pos2var]; pos2var++) { //for each position of the input barcode
        for (i = 0; i < brcd_len; i++) {     //iterate through barcode sequence
            if (i == pos2var) {              //if index i is being varied in this interation
                switch (src_ctrg->csq[i]) {  //substitute src_ctrg->csq[i] with the other 3 bases in var array
                    case 'A':
                        var[0][i] = 'T';
                        var[1][i] = 'G';
                        var[2][i] = 'C';
                        break;
                    case 'T':
                        var[0][i] = 'A';
                        var[1][i] = 'G';
                        var[2][i] = 'C';
                        break;
                    case 'G':
                        var[0][i] = 'A';
                        var[1][i] = 'T';
                        var[2][i] = 'C';
                        break;
                    case 'C':
                        var[0][i] = 'A';
                        var[1][i] = 'T';
                        var[2][i] = 'G';
                        break;
                    default:
                        printf("printSubs: error - unexpected character %c in barcode sequence. aborting...\n", src_ctrg->csq[i]);
                        abort();
                        break;
                };
            } else { //not varying index i, copy native base to var array
                var[0][i] = src_ctrg->csq[i];
                var[1][i] = src_ctrg->csq[i];
                var[2][i] = src_ctrg->csq[i];
            }
        }
        
        //add terminating null characters to each var array
        var[0][i] = '\0';
        var[1][i] = '\0';
        var[2][i] = '\0';
        
        //check that barcode was read completely
        if (i != brcd_len) {
            printf("mk_SUB_trgts: error - unexpected barcode length. aborting...\n");
            abort();
        }
        
        //store substitution variants as binary encoded sequences in compact targets
        for (i = 0; i < 3; i++) {
            set_BC_val(&ctrg[*ctrg_cnt], &BC_val[*ctrg_cnt], NULL, src_ctrg, MUTANT_BRCD); //set ptr to BC_val structure
            ctrg[*ctrg_cnt].bid = src_ctrg->bid | (*mutCode)++; //generate and store barcode id
            seq2bin_long(var[i], &ctrg[*ctrg_cnt].bsq, 1);      //store binary-encoded barcode sequence
            BC_val[*ctrg_cnt].typ = SUB;                        //set barcode type to substitution
            
            //allocate memory for storing character-encoded barcode
            if ((ctrg[*ctrg_cnt].csq = malloc((strlen(var[i])+1) * sizeof(*ctrg[*ctrg_cnt].csq))) == NULL) {
                printf("mk_SUB_trgts: error - memory allocation for character-encoded target sequence failed. aborting...\n");
                abort();
            }
            strcpy(ctrg[*ctrg_cnt].csq, var[i]); //store character-encoded barcode
            
            (*ctrg_cnt)++; //increment compact target count
        }
        sub_trgts_made += 3; //increment number of substitution targets by three
    }
    
    if (((*mutCode)-1) != MAX_SUB_CODE) {
        printf("mk_SUB_trgts: error - incorrect maximum substitution code (%llu). expected %d. aborting...\n", (long long unsigned int)((*mutCode)-1), MAX_SUB_CODE);
        abort();
    }
    
    return sub_trgts_made; //return number of substitution targets that were made
}

/* mk_INS_trgts: generate single insertion targets for input barcode */
int mk_INS_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode)
{
    int ins_trgts_made = 0;  //number of insertion targets that were made
    
    int pos2var = 0; //current position in target to add insertion
    int i = 0;       //general purpose index
    int j = 0;       //general purpose index
    
    char var[4][MAX_LINE+1] = {{0}}; //array for constructing insertion targets
    
    if (*mutCode != MIN_INS_CODE) {
        printf("mk_INS_trgts: error - incorrect minimum insertion code (%llu). expected %d. aborting...\n", (long long unsigned int)(*mutCode), MIN_INS_CODE);
        abort();
    }
    
    for (pos2var = 0; pos2var < brcd_len && src_ctrg->csq[pos2var]; pos2var++) { //for each position of the input barcode
        
        //iterate through barcode sequence
        //j starts at 1 because first input barcode character is always removed
        for (i = 0, j = 1; i < brcd_len && j <= brcd_len; i++) {
            if (i == pos2var) {  //if index i is being varied in this iteration, insert nucleotide at this position
                var[0][i] = 'A';
                var[1][i] = 'T';
                var[2][i] = 'G';
                var[3][i] = 'C';
            } else {             //otherwise, copy input sequence
                var[0][i] = src_ctrg->csq[j];
                var[1][i] = src_ctrg->csq[j];
                var[2][i] = src_ctrg->csq[j];
                var[3][i] = src_ctrg->csq[j];
                j++;
            }
        }
        
        //add terminating null character
        var[0][i] = '\0';
        var[1][i] = '\0';
        var[2][i] = '\0';
        var[3][i] = '\0';
                
        //check that input and variant barcode indices have reached full length
        if (i != brcd_len || j != brcd_len) {
            printf("mk_INS_trgts: error - unexpected barcode length. aborting...\n");
            abort();
        }
        
        //store insertion variants as binary encoded sequences in compact targets
        for (i = 0; i < 4; i++) {
            set_BC_val(&ctrg[*ctrg_cnt], &BC_val[*ctrg_cnt], NULL, src_ctrg, MUTANT_BRCD); //set ptr to BC_val structure
            ctrg[*ctrg_cnt].bid = src_ctrg->bid | (*mutCode)++;  //generate and store barcode id
            seq2bin_long(var[i], &ctrg[*ctrg_cnt].bsq, 1);       //store binary-encoded barcode sequence
            BC_val[*ctrg_cnt].typ = INS;                         //set barcode type to insertion
            
            //allocate memory for storing character-encoded barcode
            if ((ctrg[*ctrg_cnt].csq = malloc((strlen(var[i])+1) * sizeof(*ctrg[*ctrg_cnt].csq))) == NULL) {
                printf("mk_INS_trgts: error - memory allocation for character-encoded target sequence failed. aborting...\n");
                abort();
            }
            strcpy(ctrg[*ctrg_cnt].csq, var[i]); //store character-encoded barcode
            
            (*ctrg_cnt)++; //increment compact target count
        }
        ins_trgts_made += 4; //incrmement number of insertion targets by four
    }
    
    if (((*mutCode)-1) != MAX_INS_CODE) {
        printf("mk_INS_trgts: error - incorrect maximum insertion code (%llu). expected %d. aborting...\n", (long long unsigned int)((*mutCode)-1), MAX_INS_CODE);
        abort();
    }
    
    return ins_trgts_made; //return number of insertion targets made
}

/* mk_DEL_trgts: generate single deletion targets for input barcode */
int mk_DEL_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode, char upstrm_nt)
{
    int del_trgts_made = 0;  //number of deletion targets made
    
    int pos2var = 0; //current position in target to delete
    int i = 0;       //general purpose index
    int j = 0;       //general purpose index
    
    char var[MAX_LINE+1]; //array for constructing deletion targets
    
    if (*mutCode != MIN_DEL_CODE) {
        printf("mk_DEL_trgts: error - incorrect minimum deletion code (%llu). expected %d. aborting...\n", (long long unsigned int)(*mutCode), MIN_DEL_CODE);
        abort();
    }
    
    for (pos2var = 0; pos2var < brcd_len && src_ctrg->csq[pos2var]; pos2var++) { //for each position of the input barcode
        
        i = 0;                //set i to zero
        var[i++] = upstrm_nt; //set first variant base to the nt immediately upstream, since one nt will be deleted
        
        for (j = 0; i <= brcd_len && j < brcd_len; j++) { //iterate through barcode sequence
            if (j == pos2var) {                           //if pos2var is reached, do nothing so the nt is deleted
                ;
            } else {                       //othewise,
                var[i] = src_ctrg->csq[j]; //copy input barcode nt
                i++;
            }
        }
        var[i] = '\0'; //add terminating null
        
        //check that input and variant barcode indices have reached full length
        if (i != brcd_len || j != brcd_len) {
            printf("mk_DEL_trgts: error - unexpected barcode length. aborting...\n");
            abort();
        }
        
        //store deletion variant as binary encoded sequences in compact target
        set_BC_val(&ctrg[*ctrg_cnt], &BC_val[*ctrg_cnt], NULL, src_ctrg, MUTANT_BRCD); //set ptr to BC_val structure
        ctrg[*ctrg_cnt].bid = src_ctrg->bid | (*mutCode)++;  //generate and store barcode id
        seq2bin_long(var, &ctrg[*ctrg_cnt].bsq, 1);          //store binary-encoded barcode sequence
        BC_val[*ctrg_cnt].typ = DEL;                         //set barcode type to deletion
        
        //allocate memory for storing character-encoded barcode
        if ((ctrg[*ctrg_cnt].csq = malloc((strlen(var)+1) * sizeof(*ctrg[*ctrg_cnt].csq))) == NULL) {
            printf("mk_DEL_trgts: error - memory allocation for character-encoded target sequence failed. aborting...\n");
            abort();
        }
        strcpy(ctrg[*ctrg_cnt].csq, var); //store character-encoded barcode
        
        (*ctrg_cnt)++;     //increment compact target count
        del_trgts_made++;  //increment number of deletion targets made by one
    }
    
    if (((*mutCode)-1) != MAX_DEL_CODE) {
        printf("mk_DEL_trgts: error - incorrect maximum deletion code (%llu). expected %d. aborting...\n", (long long unsigned int)((*mutCode)-1), MAX_DEL_CODE);
        abort();
    }
    
    return del_trgts_made; //return number of deletion targets made
}

/* get_target_type: determine target type using mutcode */
char * get_target_type(uint64_t mutcode, int * type_val)
{
    static char nat_type[4] = {"NAT"};
    static char sub_type[4] = {"SUB"};
    static char ins_type[4] = {"INS"};
    static char del_type[4] = {"DEL"};
    
    if (!mutcode) {
        *type_val = NAT;
        return nat_type;
        
    } else if (mutcode >= MIN_SUB_CODE && mutcode <= MAX_SUB_CODE) {
        *type_val = SUB;
        return sub_type;
        
    } else if (mutcode >= MIN_INS_CODE && mutcode <= MAX_INS_CODE) {
        *type_val = INS;
        return ins_type;
        
    } else if (mutcode >= MIN_DEL_CODE && mutcode <= MAX_DEL_CODE) {
        *type_val = DEL;
        return del_type;
        
    } else {
        printf("error: mutcode value is out of expected bounds. aborting...\n");
        abort();
    }
}

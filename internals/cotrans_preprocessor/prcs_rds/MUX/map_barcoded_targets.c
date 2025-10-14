//
//  map_barcoded_targets.c
//  
//
//  Created by Eric Strobel on 10/14/25.
//

#include <stdio.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../../utils/io_management.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "../../MUX_trgt_gen/mk_MUX_trgts.h"
#include "../../MUX_trgt_gen/set_barcoded_compact_target.h"


#include "map_barcoded_targets.h"

/* mk_htbl_MUX: makes compact target hash table */
/* hash table has linked list buckets for possible collisions */
void mk_htbl_MUX(compact_h_node ** htbl_MUX, compact_h_node_bank * bank, compact_target * ctrg, int count, mapping_metrics * met)
{
    extern uint64_t mutcode_mask; //mask to used to isolate barcode mutation code
    
    compact_h_node **p_rdnd = NULL; //pointer for h_node handling
    
    int i = 0;              //general purpose index
    int j = 0;              //general purpose index
    int new_node = 0;       //counts number of nodes that are assigned a target
    int redundant = 0;      //counts number of redundant targets (same seq as prev target)
    int collisions = 0;     //counts number of collisions between barcodes
    int ntv_collisions = 0; //counts number of native barcode collisions
    int blacklisted = 0;    //counts number of blacklisted mutant barcodes
    
    int error = 0; //flag that error occurred and program should be aborted
    
    uint64_t hash = 0;  //hash key
    
    char csq1[MAX_LINE+1] = {0};  //array for storing character-encoded barcode sequence
    char csq2[MAX_LINE+1] = {0};  //array for storing character-encoded barcode sequence
    
    opt_BC * p_BC_val_1 = NULL; //pointer to barcode target optional values
    opt_BC * p_BC_val_2 = NULL; //pointer to barcode target optional values
    
    char * type_str1 = NULL; //pointer to barcode type string for sequence 1
    char * type_str2 = NULL; //pointer to barcode type string for sequence 2
    int type_val1 = -1;      //sequence 1 barcode type value
    int type_val2 = -1;      //sequence 2 barcode type value
    
    for (i = 0; i < count; i++) {  //perform loop for each 3' end target
        
        hash = hash_brcd_trgt(&ctrg[i].bsq);                      //hash binary encoded sequence
        p_rdnd = srch_ctrg_htbl(&ctrg[i].bsq, hash, htbl_MUX, 0); //search hash table
        
        if ((*p_rdnd) == NULL) {                    //no existing hash table node for target sequence
            (*p_rdnd) = &bank->chn[bank->count++];  //assign node from hash node bank
            if (bank->count == BLOCK_SIZE) {        //check that bank was not filled
                extend_ch_bank(bank);               //extend bank if needed
            }
            (*p_rdnd)->ctrg = &(ctrg[i]);           //set node to point to current target
            new_node++;                             //increment new_node counter
            
        } else {
            ctrg[i].mul = 1; //set flag that current target is redudant with previous target
            redundant++;     //increment redundant counter
            
            //check whether redundant barcode mutants are linked to different reference barcodes
            if (check_brcd_diff((*p_rdnd)->ctrg, &ctrg[i])) {     //if brcds map to different ref brcds
                collisions++;                                     //increment collision counter
                
                bin2seq(csq1, &(*p_rdnd)->ctrg->bsq, MAX_LINE+1); //convert barcode 1 to string
                bin2seq(csq2, &ctrg[i].bsq, MAX_LINE+1);          //convert barcode 2 to string
                
                p_BC_val_1 = (opt_BC *)(*p_rdnd)->ctrg->opt; //set pointer to barcode 1 optional values
                p_BC_val_2 = (opt_BC *)ctrg[i].opt;          //set pointer to barcode 2 optional values
                
                type_str1 = get_target_type((*p_rdnd)->ctrg->bid & mutcode_mask, &type_val1);
                type_str2 = get_target_type(ctrg[i].bid & mutcode_mask, &type_val2);
                
                //report collision
                printf("\nCOLLISION:\n");
                
                if (!strcmp(csq1, csq2)) { //collided barcode sequences match (sanity check)
                    //printf("id1:  ");
                    //printbin((*p_rdnd)->ctrg->bid);
                    //printf("id2:  ");
                    //printbin(ctrg[i].bid);
                    printf("ref1: %s\nref2: %s\n", p_BC_val_1->ref->csq, p_BC_val_2->ref->csq);
                    printf("trg1: %s (%s)\ntrg2: %s (%s)\n\n", csq1, type_str1, csq2, type_str2);
                    
                } else { //collided barcode sequence do not match - this shouldn't be possible
                    printf("mk_htbl_mux: error - collision of non-identical barcodes. this should not be possible.\n");
                    printf("brcd: %s\n", csq1);
                    printf("brcd: %s\n", csq2);
                    printf("aborting...\n");
                    abort();
                }
                
                
                //TODO: since native barcodes are first and all native barcodes are guaranteed to be sufficiently distant from each other, I think only (*p_rdnd)->ctrg->bid needs to be tested here. Is there a case where ctrg[i].bid could be a native barcode?
                //test whether collision involved reference barcode
                if (!((*p_rdnd)->ctrg->bid & mutcode_mask) || !(ctrg[i].bid & mutcode_mask)) { //if either brcd was NAT

                    ntv_collisions++; //increment number of native barcode collisions
                    
                    //it should not be possible for a collision in which the second barcode was native to occur.
                    //native barcodes are added to the hash table first, and all native barcodes should be distinct
                    if (!(ctrg[i].bid & mutcode_mask)) { //if second barcode was native, throw error and abort
                        printf("\nmk_htbl_MUX: error - native collision in which second target was native barcode. this should not be possible. aborting after hash table generation...\n");
                        error = 1;
                    } else {                //otherwise
                        ctrg[i].bl++;       //blacklist second barcode
                        met->blacklisted++; //increment number of blacklisted barcodes
                        //note mul does not need to be set here since it is set above for all redundant
                        //targets regardless of whether a collision occurred
                    }
                    
                } else { //both barcode targets were mutants, set prev target as redundant, blacklist both
                    (*p_rdnd)->ctrg->mul = 1; //set previous target as redundant
                    (*p_rdnd)->ctrg->bl++;    //blacklist previous target
                    ctrg[i].bl++;             //blacklist current target
                    //as above, mul does not need to be set here since it is set above for all redundant
                    //targets regardless of whether a collision occurred
                }
            }
        }
    }
    
    //print hash table generation details
    printf("\n********************  hash table generation  ********************\n");
    printf("%6d barcode targets were assessed\n", i);
    printf("%6d barcode target sequences were assigned a node\n", new_node);
    printf("%6d redundant barcode target sequences were not assigned a node\n\n\n", redundant);
    
    if (error) { //if a serious error occurred above, throw error and abort
        printf("aborting due to fatal error. see message above for cause.\n");
        abort();
    }
    
    //prompt user regarding whether analysis should proceed
    //TODO: can likely remove this. if a collision occurs with a native barcode the program is stopped. collisions between mutant barcodes are unlikely to impact the results of the analysis since reads that map to mutant barcodes should be rare. instead, consider generating a file that reports all collisions so that there is a record
    char usr_resp[4] = {0};         //storage for user response
    char discard[MAX_LINE+1] = {0}; //array for flushing stdin
    int resp_provided = 0;          //flag that valid reponse was provided
    
    if (collisions) {
        printf("%d collision(s) between native/substitution/insertion/deletion barcode sequences were detected. proceed?\n", collisions);
        
        while (!resp_provided) {
            scanf("%3s", usr_resp);
            
            if (!strcmp(usr_resp, "yes")) {
                printf("proceeding with variant demultiplexing\n\n");
                resp_provided = 1;
                
            } else if (!strcmp(usr_resp, "no")) {
                printf("aborting variant demultiplexing\n\n");
                abort();
                
            } else {
                printf("invalid response. please enter \"yes\" or \"no\".\n");
                get_line(discard, stdin);
                
            }
        }
    }
}

/* hash_brcd_trgt: generates hash key for binary encoded sequence (up to 32 nt/64 bits */
uint64_t hash_brcd_trgt(binary_seq * bsq)
{
    return *bsq->sq % TABLE_SIZE;
}

/* check_brcd_diff: checks whether redundant barcode targets are linked to different barcode variants*/
int check_brcd_diff(compact_target * old, compact_target * new)
{
    if ((new->bid >> MUTCODE_BITS) != (old->bid >> MUTCODE_BITS)) {
        return 1;
    } else {
        return 0;
    }
}
/* map_brcd: map barcode to target using hash table */
compact_target * map_brcd(char * brcd_str, char * rc_brcd_str, compact_h_node **htbl_MUX, compact_target ** mpd_trg, mapping_metrics * met)
{
    extern struct testdata_MUX_vars testdata_MUX; //structure containing test data read analysis variables
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    compact_h_node ** p_rdnd = NULL; //pointer to hash table node pointer
    compact_target * ref_trg = NULL; //pointer to reference compact target structure
    
    opt_BC * BC_val  = NULL; //pointer to current barcode optional target value structure
    opt_BC * ref_val = NULL; //pointer to reference barcode optional target value structure
    
    binary_seq bsq = {0}; //variable for storing binary-encoded barcode sequence
    if ((bsq.sq = calloc(1, sizeof(*(bsq.sq)))) == NULL) {
        printf("map_brcd: error - failed to allocate memory for binary-encoded barcode sequence. aborting\n");
        abort();
    }
    
    uint64_t hash = 0; //hash value
    
    seq2bin_long(rc_brcd_str, &bsq, 1); //generate binary-encoded barcode sequence
    hash = hash_brcd_trgt(&bsq);        //hash binary-encoded barcode sequence
    
    //TODO: for now, trace search is off. add as debug option later
    p_rdnd = srch_ctrg_htbl(&bsq, hash, htbl_MUX, 0); //search hash table for match
    
    if ((*p_rdnd) != NULL) {                              //if found barcode match
        met->hits++;                                      //increment number of hits
        if (!strcmp(rc_brcd_str, (*p_rdnd)->ctrg->csq)) { //sanity check that barcode sequences match
            met->matches++;                               //if match, increment number of barcode sequence matches
        } else {                                          //otherwise, throw error and abort
            printf("map_brcd: error query string does not match has table entry. aborting...\n");
            abort();
        }
        
        if ((*p_rdnd)->ctrg->bl) { //if target is blacklisted,
            met->bl_mapped++;      //track mapped blacklisted read count
            return NULL;           //then return NULL
            
        } else { //otherwise track mapped read and return pointer to reference target
            met->mapped++;                           //track mapped read
            *mpd_trg = (*p_rdnd)->ctrg;              //set mapped target pointer to the current node's target
            
            BC_val = (opt_BC *)(*p_rdnd)->ctrg->opt; //set pointer to current target optional values
            BC_val->cnt++;                           //increment number of reads that mapped to current target
            switch (BC_val->typ) {                   //track number of NAT/SUB/INS/DEL reads
                case NAT: met->nat_cnt++; break;
                case SUB: met->sub_cnt++; break;
                case INS: met->ins_cnt++; break;
                case DEL: met->del_cnt++; break;
                default:
                    printf("map_brcd: unrecognized barcode mutant type. aborting...\n");
                    abort();
                    break;
            }
        
            ref_trg = BC_val->ref;            //set pointer to reference target of current mapped target
            ref_val = (opt_BC *)ref_trg->opt; //set pointer to reference target optional values
            ref_val->mpd++;                   //increment number of reads that mapped to the reference target
            free(bsq.sq);                     //free memory
            return ref_trg; //return reference target pointer
        }
    } else {             //otherwise
        met->unmapped++; //increment number of unmapped reads
        free(bsq.sq);    //free memory
        return NULL;     //return NULL
    }
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

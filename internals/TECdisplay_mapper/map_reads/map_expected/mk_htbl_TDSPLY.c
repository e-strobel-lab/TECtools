//
//  mk_htbl_TDSPLY.c
//  
//
//  Created by Eric Strobel on 7/13/26.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

#include "../../../utils/debug.h"

#include "../../../seq_utils/seq2bin_hash.h"

#include "./get_key.h"


#include "mk_htbl_TDSPLY.h"


/* mk_htbl_TDSPLY: construct hash table from target structure */
int mk_htbl_TDSPLY(h_node **htbl, h_node_bank *bank, target *trgts, target *refs, target_params *trg_prms)
{
    extern int debug;
    
    h_node **p_refnd = NULL; //pointer for h_node handling from reference key look up
    h_node **p_altnd = NULL; //pointer for h_node handling during alternate key search
    h_node  *null_nd = NULL; //null h_node pointer for initializing p_altnd in loop below
    
    int i = 0; //general purpose index
    int r = 0; //reference target index
    
    char altkey[MAX_LINE] = {0}; //array for storing alternate reference target keys
    int fnd_altMatch = 0;        //flag that alternate reference target match was found
    
    int new_node = 0;  //number of nodes that are assigned a target
    int redundant = 0; //number of redundant targets (same seq as prev target)
    
    /* targets files can be generated using multiple reference sequences. to identify
     redundant targets from different reference sequences, hash table searches are
     performed for each target using a key generated from its reference sequence and
     keys generated from all other reference sequences that were present in the targets
     file. if the target is non-redundant, it is assigned a hash table node using the
     key that was generated from its reference sequence. redundant targets are ignored.
     */
    
    for (i = 0; i < trg_prms->t_cnt; i++) {      //perform loop for each target
        if (debug) {printf("target %5d: %s\nsource reference key:\n", i, trgts[i].id);}
        p_refnd = srch_htbl(trgts[i].key, htbl); //search hash table with key from source reference sequence
        
        //search hash table with keys generated using other
        //reference sequences until all keys have been tested
        //or a match is found. the source reference sequences
        //is tested again here, but this doesn't impact
        //performance in a meaningful way
        if (debug) {printf("all reference key check:\n");}
        for (r = 0, p_altnd = &null_nd, fnd_altMatch = 0; r < trg_prms->r_cnt && !fnd_altMatch; r++) {
            p_altnd = &null_nd;                                             //set p_altnd to NULL
            get_key(altkey, trgts[i].sq, NULL, NULL, &refs[r], TARGET_KEY, MAX32CHAR_TARGET); //generate key using alt reference seq
            
            if ((p_altnd = srch_htbl(altkey, htbl))) {                //search hash table with alternate reference key
                if ((*p_altnd) != NULL) {                             //if a node match is found
                    if (!strcmp(trgts[i].sq, (*p_altnd)->trg->sq)) {  //check for full sequence match
                        fnd_altMatch = 1;                             //set flag that full seuence match was found
                    }
                }
            }
            if (debug) {printf("ref%d\t%s key=%s\n", r, (*p_altnd == NULL) ? "NULL " : "MATCH", altkey);}
        }
        if (debug) {printf("\n");}
        
        /* TODO: need to set up to handle targets with identical hash keys that are distinct */
        
        if ((*p_refnd) == NULL && !fnd_altMatch) {      //no existing hash table node for target sequence
            (*p_refnd) = &bank->hn[bank->count++];      //assign node from hash node bank
            (*p_refnd)->trg = &(trgts[i]);              //set node to point to current target
            
            if (bank->count == BLOCK_SIZE) {            //check whether bank was filled
                extend_h_bank(bank);                    //if filled, extend bank
            }
            
            new_node++;                                 //increment new_node counter
        } else {
            trgts[i].mul = 1; //set flag that current target is redudant with previous target
            redundant++;      //increment redundant counter
        }
    }
    printf("\n********************  hash table generation  ********************\n");
    printf("%6d targets were assessed\n", i);
    printf("%6d target sequences were assigned a node\n", new_node);
    printf("%6d redundant target sequences were not assigned a node\n\n\n", redundant);
    
    trg_prms->nr_cnt = new_node;  //set non-redundant node count
    trg_prms->rdndnt = redundant; //set redundant node count
    
    return new_node; //return number of non-redundant targets
}

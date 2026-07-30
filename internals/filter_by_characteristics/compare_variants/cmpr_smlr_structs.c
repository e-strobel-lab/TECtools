//
//  cmpr_smlr_structs.c
//  
//
//  Created by Eric Strobel on 6/29/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../set_seq_attributes/parse_descriptor_file.h"

#include "./mk_dot_bracket_htbl.h"

#include "cmpr_smlr_structs.h"

/* cmpr_smlr_structs: perform a systematic comparison of variants with closely-related structures to identify variants that exhibit large effects on function */
void cmpr_smlr_structs(sequence_attributes * sq_att, descriptor * des, int seq_cnt, int des_cnt, char * des_nm)
{
    /*
     a hash_table for looking up variants that exhibit a particular secondary structure is generated as follows:
     
     1. in cmpr_smlr_structs:
        - allocate compact_target, opt_db, and sp_list pointer memory for every
          structProps associated with the descriptor being used for comparisons

     2. in set_db_compact_target:
        - generate a compact_target for every structProps structure
        - associate each structProps structure with its corresponding compact_target
        - add each structProps structure to sp_list

     3. in mk_dot_bracket_htbl:
        - map every compact_target (generated from structProps) to the hash table.
          -the number of subsequence compact_targets that match the first instance of a structure are using mul
           (in the first instance ctrg) and blacklisted so that only the first instance of a structure is
           considered during downstream processing
     
        - using init_opt_db, structProps memory associated with an opt_db struct is allocated for every
          structProps structure that matches the first instance of a structure.
     
        - using map_structProps_2_htbl, each structProps is associated with its corresponding compact_target
          such that all variants that exhibit a given structures are associated by the compact target
     */
    
    
    compact_target *ctrg = NULL;    //pointer for allocating compact_target structs
    opt_db * db_val = NULL;         //pointer for allocating opt_db structs
    
    structProps ** sp_list = NULL;  //pointer for generating structProps structure list
    
    int ctrg_cnt = 0;         //number of compact targets
    
    int i = 0;                //general purpose index
    int fnd_des = 0;          //flag that descriptor to be assessed was found
    int des_ix = 0;           //index of descriptor to be assessed
    int des_typ = TYPE_INIT;  //type of descriptor to be assessed
    
    //search array of descriptors to identify match to des_nm, which is the descriptor to be assessed
    for (i = 0, fnd_des = 0, des_ix = 0; i < des_cnt && !fnd_des; i++) {
        if (!strcmp(des[i].nm, des_nm)) { //if a match to des_nm is found
            ctrg_cnt = des[i].tot_sp;     //set compact target count to the total number of structProps associated with the descriptor
            des_ix = i;                   //set des_ix to index of the current descriptor match
            des_typ = des[i].typ;         //set des_typ to the type of the current descriptor
            fnd_des = 1;                  //set flag that match to current descriptor was found
        }
    }
    
    if (!fnd_des) { //if no match to des_nm was found, throw error and abort
        printf("cmpr_smlr_structs: error - could not find descriptor with name '%s'. aborting...\n", des_nm);
        abort();
    }
    
    printf("cmpr_smlr_structs: tot_sp=%d\n", ctrg_cnt);
    
    //allocate memory for compact targets
    if ((ctrg = calloc(ctrg_cnt, sizeof(*ctrg))) == NULL) {
        printf("cmpr_smlr_structs: error - target memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //allocate memory for optional values
    if ((db_val = calloc(ctrg_cnt, sizeof(*db_val))) == NULL) {
        printf("cmpr_smlr_structs: error - opt values memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //allocate memory for structProps list
    if ((sp_list = calloc(ctrg_cnt, sizeof(*sp_list))) == NULL) {
        printf("cmpr_smlr_structs: error - structProps list memory allocation failed\n");
        abort(); //TODO: abort or return error code?
    }
    
    //generate compact_targets from structProps structures
    set_db_compact_target(ctrg, db_val, sp_list, ctrg_cnt, sq_att, seq_cnt, des_ix, des_typ);
    
    /********* hash table initialization and construction **********/
    compact_h_node **htbl = NULL;                  //hash table root
    compact_h_node_bank hn_bank = {NULL, NULL, 0}; //bank for hash table nodes
    
    hn_bank.count = 0; //initialize hash table node count to zero
    
    //allocate TABLE_SIZE hash table node pointers
    if ((htbl = calloc(TABLE_SIZE, sizeof(*htbl))) == NULL) {
        printf("cmpr_smlr_structs: error - hash table memory allocation failed. aborting...\n");
        abort();
    }
    
    //allocate BLOCK_SIZE hash table nodes
    if ((hn_bank.chn = calloc(BLOCK_SIZE, sizeof(*(hn_bank.chn)))) == NULL) {
        printf("cmpr_smlr_structs: error - hash table node memory allocation failed. aborting...\n");
        abort();
    }
    
    //generate htbl using dot-bracket structures
    compact_h_node_bank *crrnt_hn_bank = &hn_bank;    //pointer for handling hash node bank
    mk_dot_bracket_htbl(htbl, crrnt_hn_bank, ctrg, ctrg_cnt, sp_list/*, target_params * trg_prms*/);
    
    return;
}

/* set_db_compact_target: generate compact_targets from structProps structures and associate the structProps struct with the compact_target */
void set_db_compact_target(compact_target * ctrg, opt_db * db_val, structProps ** sp_list, int ctrg_cnt, sequence_attributes * sq_att, int seq_cnt, int des_ix, int des_typ)
{
    int i = 0;      //general purpose index
    int j = 0;      //general purpose index
    int c = 0;      //compact_target index
    int sp_cnt = 0; //number of structProps associated with the current attribute
    
    proximal_deltaG * p_prx_dG = NULL; //pointer for dereferencing proximal_deltaG attributes
    distal_deltaG * p_dst_dG = NULL;   //pointer for dereferencing distal_deltaG attributes
    structProps * p_sp = NULL;         //pointer to current structProps structure
        
    for (i = 0, c = 0; i < seq_cnt; i++) { //for every variant sequence
        
        //set values for structProps of current variant sequence
        switch (des_typ) {
            case PRX_DG:                                             //descriptor is PRX_DG
                p_prx_dG = (proximal_deltaG *)sq_att[i].att[des_ix]; //dereference attribute as PRX_DG
                p_sp = &p_prx_dG->sp;                                //set pointer to root structProps
                sp_cnt = p_prx_dG->sp_cnt;                           //set number of structProps associated with current attribute
                break;
                
            case DST_DG:                                             //descriptor is DST_DG
                p_dst_dG = (distal_deltaG *)sq_att[i].att[des_ix];   //dereference attribute as DST_DG
                p_sp = &p_dst_dG->sp;                                //set pointer to root structProps
                sp_cnt = p_dst_dG->sp_cnt;                           //set number of structProps associated with current attribute
                break;
                
            default: //unrecognized descriptor type, throw error and abort
                printf("set_db_compact_target: error - unrecognized descriptor type. aborting...\n");
                abort();
                break;
        }
    
        //generate compact_target from current structProps and
        //and associate current structProps with the compact_target
        for (j = 0; j < sp_cnt; j++) {    //for every structProps associated with the current attribute
                        
            if (j > 0) {                  //if not processing first structProps
                if (p_sp->nxt != NULL) {  //if next pointer to next structProps in linked list is not null
                    p_sp = p_sp->nxt;     //set p_sp to next structProps struct
                    
                } else { //fewer than sp_cnt structProps in linked list, throw error and abort
                    printf("set_db_compact_target: error - fewer than expected structProps structures in the linked list. aborting...\n");
                    abort();
                }
            }
            
            //generate compact_target using dot-bracket structure from current structProps
            if (c < ctrg_cnt) { //check that current ctrg is not exceed allocated ctrgs
                
                //allocate memory for storing dot-bracket sequence in compact_target structure
                if ((ctrg[c].csq = malloc((strlen(p_sp->db)+1) * sizeof(*ctrg[c].csq))) == NULL) {
                    printf("set_db_compact_target: error - failed to allocate compact target char seq memory. aborting...\n");
                    abort();
                }
                
                strcpy(ctrg[c].csq, p_sp->db);                                  //store dot-bracket structure string in compact_target
                dot_bracket_to_bin_hash(ctrg[c].csq, &ctrg[c].bsq, MAX_BLOCKS); //convert dot-bracket structure to binary
                ctrg[c].opt = &db_val[c];                                       //point opt to current opt_db struct
                sp_list[c] = p_sp;                                              //add current structProps to sp_list
                c++;                                                            //increment ctrg index
                
            } else { //not enough compact_targets, throw error and abort
                printf("set_db_compact_target: error - not enough compact targets to store all structure properties entries. aborting...\n");
                abort();
            }
        }
    }
    
    return;
}

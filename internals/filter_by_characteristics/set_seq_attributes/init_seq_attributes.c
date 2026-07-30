//
//  init_seq_attributes.c
//  
//
//  Created by Eric Strobel on 11/21/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

#include "init_seq_attributes.h"

/* init_seq_attributes: initialize sequence attributes structure */
void init_seq_attributes(sequence_attributes ** sq_att, descriptor * des, descriptor_bank * des_bnk, int seq_cnt, int des_cnt, int * typ_cnt)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    if (typ_cnt[NUC_ID]) { //if NUC_ID types were detected, allocate NUC_ID attribute bank memory
        if ((des_bnk->nuc_id = calloc(typ_cnt[NUC_ID] * seq_cnt, sizeof(*des_bnk->nuc_id))) == NULL) {
            printf("init_seq_attributes: error - failed to allocate attribute memory. aborting...\n");
            abort();
        }
    }
    
    if (typ_cnt[PRX_DG]) { //if PRX_DG descriptor types were detected, allocate PRX_DG attribute bank memory
        if ((des_bnk->prx_dG = calloc(typ_cnt[PRX_DG] * seq_cnt, sizeof(*des_bnk->prx_dG))) == NULL) {
            printf("init_seq_attributes: error - failed to allocate attribute memory. aborting...\n");
            abort();
        }
    }
    
    if (typ_cnt[DST_DG]) { //if DST_DG descriptor types were detected, allocate DST_DG attribute bank memory
        if ((des_bnk->dst_dG = calloc(typ_cnt[DST_DG] * seq_cnt, sizeof(*des_bnk->dst_dG))) == NULL) {
            printf("init_seq_attributes: error - failed to allocate attribute memory. aborting...\n");
            abort();
        }
    }
    
    if (typ_cnt[SS_LEN]) { //if SS_LEN descriptor types were detected, allocate SS_LEN attribute bank memory
        if ((des_bnk->ss_len = calloc(typ_cnt[SS_LEN] * seq_cnt, sizeof(*des_bnk->ss_len))) == NULL) {
            printf("init_seq_attributes: error - failed to allocate attribute memory. aborting...\n");
            abort();
        }
    }
        
    //allocate and assign attribute memory for all input sequences
    for (i = 0; i < seq_cnt; i++) {
        
        //allocate attribute memory for each each descriptor
        if (((*sq_att)[i].att = calloc(des_cnt, sizeof(*((*sq_att)[i].att)))) == NULL) {
            printf("init_seq_attributes: error - failed to allocate attribute memory. aborting...\n");
            abort();
        }
        
        (*sq_att)[i].des = des; //set pointer to descriptor structure
        
        //assign each attribute pointer an attribute structure from the corresponding bank
        for (j = 0; j < des_cnt; j++) {
            switch (des[j].typ) {
                case NUC_ID: (*sq_att)[i].att[j] = &des_bnk->nuc_id[des_bnk->cnt[NUC_ID]++]; break;
                case PRX_DG: (*sq_att)[i].att[j] = &des_bnk->prx_dG[des_bnk->cnt[PRX_DG]++]; break;
                case DST_DG: (*sq_att)[i].att[j] = &des_bnk->dst_dG[des_bnk->cnt[DST_DG]++]; break;
                case SS_LEN: (*sq_att)[i].att[j] = &des_bnk->ss_len[des_bnk->cnt[SS_LEN]++]; break;
                    
                default:
                    printf("init_seq_attributes: error - unrecognized descriptor type. aborting...\n");
                    abort();
                    break;
            }
        }
    }
    
    //check that the number of assigned attribute structures matches expectations
    for (i = 0; i < DSCRPTR_TYPES; i++) {
        if (des_bnk->cnt[i] != (typ_cnt[i] * seq_cnt)) {
            printf("init_seq_attributes: error - number of assigned attribute structures does not match calculated number. aborting...\n");
            abort();
        }
    }
    
    return;
}

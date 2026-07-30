//
//  set_compact_target.c
//  
//
//  Created by Eric Strobel on 7/6/26.
//

#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"
#include "../../../seq_utils/basemap.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

#include "./parse_vmt_trgts.h"

#include "../../../filter_by_characteristics/filter_by_characteristics_defs.h"
#include "../../../filter_by_characteristics/filter_by_characteristics_structs.h"
#include "../../../filter_by_characteristics/set_seq_attributes/parse_descriptor_file.h"

#include "set_compact_target.h"

/* set_trgt: set target values in compact target struct */
void set_compact_target(compact_target * ctrg, opt_mx_trg * trg_val, target * crnt_ref, char * trgt_id, char * trgt_sq)
{
    opt_mx_trg * p_trg_val = NULL; //pointer for dereferencing opt pointer in targets structure to opt_mx_trg
    opt_ref * p_ref_val = NULL;    //pointer for dereferencing opt pointer in targets structure to opt_ref
    
    //filter non-base chars from input seq
    char processed_seq[MAXLEN+1] = {0}; //array to store seq w/o non-base chars
    process_trgt_seq(trgt_sq, processed_seq);
    
    //allocate memory for target structure values
    if ((ctrg->cid = malloc((strlen(trgt_id)+1) * sizeof(*ctrg->cid))) == NULL) {
        printf("set_compact_target: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    if ((ctrg->csq = malloc((strlen(trgt_sq)+1) * sizeof(*ctrg->csq))) == NULL) {
        printf("set_compact_target: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    //set compact targets structure values
    strcpy(ctrg->cid, trgt_id);       //copy target name to target entry
    strcpy(ctrg->csq, processed_seq); //copy processed target sequence to target entry
    
    ctrg->opt = trg_val;                        //point opt to corresponding target values
    p_trg_val = (opt_mx_trg *)ctrg->opt;        //dereference trgts->opt to simplify code below
    p_trg_val->ref = crnt_ref;                  //set trgts->opt->ref to point to current reference target
    
    p_ref_val = (opt_ref *)p_trg_val->ref->opt;  //dereference trgts->opt->ref to simplify code below
    
    if ((p_trg_val->vb_key = calloc(p_ref_val->vb_cnt+1, sizeof(*p_trg_val->vb_key))) == NULL) {
        printf("set_compact_trgt: error - failed to allocate memory for variable base key. aborting...\n");
        abort();
    }

    get_key(p_trg_val->vb_key, ctrg->csq, NULL, NULL, p_trg_val->ref, TARGET_KEY, UNBOUNDED_TARGET); //generate key
    seq2bin_long(p_trg_val->vb_key, &ctrg->bsq, MAX_BLOCKS); //store key as binary_seq
    
    return;
}

/* set_utility_target: set standard target as utility member of compact target struct */
void set_utility_target(compact_target * ctrg, target * utl, int cnt)
{
    int i = 0; //general purpose index
    
    opt_mx_trg * p_trg_val = NULL; //pointer for dereferencing opt_mx_trg from compact_target opt member
    
    for (i = 0; i < cnt; i++) { //for every compact target
        
        ctrg[i].utl = &utl[i]; //point utl pointer to utility target
       
        p_trg_val = (opt_mx_trg *)ctrg[i].opt; //dereference target values
        
        utl[i].key = p_trg_val->vb_key; //set key pointer
        utl[i].id  = ctrg[i].cid;       //set id pointer
        utl[i].sq  = ctrg[i].csq;       //set sequence pointer
        utl[i].cnt = ctrg[i].cnt;       //set mapped read count
        utl[i].mul = ctrg[i].mul;       //set redundant flag
        utl[i].opt = ctrg[i].opt;       //set opt pointer in utl to point to opt associated with ctrg
    }
    
    return;
}

//
//  set_trgt.c
//  
//
//  Created by Eric Strobel on 9/30/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/revcomp.h"
#include "../../../seq_utils/isDNAbase.h"
#include "../../../seq_utils/isIUPACbase.h"
#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/basemap.h"

#include "./get_key.h"
#include "./parse_mx_trgts.h"

#include "set_trgt.h"

/* set_trgt: set target values in target struct */
void set_trgt(target * trgts, opt_mx_trg * trg_val, target * crnt_ref, char * trgt_id, char * trgt_sq)
{
    extern int debug;
    
    int i = 0; //general purpose index
    
    opt_mx_trg * p_trg_val = NULL;     //pointer for dereferencing opt pointer in targets structure to opt_mx_trg
    opt_ref * p_ref_val = NULL;         //pointer for dereferencing opt pointer in targets structure to opt_ref
    
    //filter non-base chars from input seq
    char processed_seq[MAXLEN+1] = {0}; //array to store seq w/o non-base chars
    if (!process_trgt_seq(trgt_sq, processed_seq)) {
        printf("set_trgt: error - non-base characters in reference %s and target %s do not match. aborting...", crnt_ref->id, trgt_id);
        abort();
    }
        
    //allocate memory for target structure values
    if ((trgts->id = malloc((strlen(trgt_id)+1) * sizeof(trgts->id))) == NULL) {
        printf("set_trgt: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    if ((trgts->sq = malloc((strlen(processed_seq)+1) * sizeof(trgts->sq))) == NULL) {
        printf("set_trgt: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    //set trgts structure values
    strcpy(trgts->id, trgt_id);       //copy target name to target entry
    strcpy(trgts->sq, processed_seq); //copy processed target sequence to target entry
    
    trgts->opt = trg_val;                        //point opt to corresponding target values
    p_trg_val = (opt_mx_trg *)trgts->opt;        //dereference trgts->opt to simplify code below
    p_trg_val->ref = crnt_ref;                   //set trgts->opt->ref to point to current reference target
    p_ref_val = (opt_ref *)p_trg_val->ref->opt;  //dereference trgts->opt->ref to simplify code below
    
    get_key(p_trg_val->vb_key, trgts->sq, NULL, NULL, p_trg_val->ref, TARGET_KEY); //generate key
    trgts->key = &(p_trg_val->vb_key[0]); //set key pointer to point to key in opt_mx_trg structure
    
    return;
}

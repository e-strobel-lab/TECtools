//
//  parse_fa_trgts.c
//  
//
//  Created by Eric Strobel on 10/7/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"

#include "./set_barcoded_compact_target.h"

#include "parse_fa_trgts.h"


/* parse_fa_trgts: stores input barcodes as targets */
void parse_fa_trgts(FILE * ifp, int trgt_ftype, void * trgts, void * trg_val, target_params * trg_prms, int data_type)
{
    uint64_t i = 0; //general purpose index
    
    opt_BC * p_opt_BC = NULL; //pointer for handling barcode target optional values
    
    char trgt_id[MAX_LINE+1] = {0};  //array for current target identifier
    char trgt_sq[MAX_LINE+1] = {0};  //current target sequence
        
    int L1 = 0; //variable for checking get fasta line 1 success
    int L2 = 0; //variable for checking get fasta line 2 success
    
    while ((L1 = get_line(trgt_id, ifp))) {
        
        //if t_cnt > xpctd and got_line was successful,
        //there are more targets than expected and not
        //enough memory was allocated.
        if (trg_prms->t_cnt > trg_prms->xpctd) {
            printf("parse_fa_trgts: error - more targets than expected in targets file. aborting...\n");
            abort();
        }
        
        L2 = get_line(trgt_sq, ifp); //get fasta line 2
        
        if (L1 == 0 || L2 == 0) { //if getting line 1 or line 2 failed, throw error and abort
            printf("parse_fa_trgts: error - reached end of barcodes file before reading expected number of barcodes. aborting...\n");
            abort();
        }
        
        set_barcoded_compact_target(&(((compact_target *)trgts)[trg_prms->t_cnt]),
                                    &(((opt_BC *)trg_val)[trg_prms->t_cnt]),
                                    NULL, trgt_id, trgt_sq, trg_prms, trgt_ftype, data_type);
        
        trg_prms->t_cnt++; //increment target count
        L1 = L2 = 0;       //zero L1 and L2
    }
    
    return;
}

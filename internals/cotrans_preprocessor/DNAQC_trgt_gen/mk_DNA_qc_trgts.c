//
//  mk_DNA_qc_trgts.c
//  
//
//  Created by Eric Strobel on 9/9/26.
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

#include "../MUX_trgt_gen/mk_MUX_trgts.h"

#include "./parse_DNA_qc_trgts.h"

#include "mk_DNA_qc_trgts.h"

/* mk_DNA_qc_trgts: generate targets for DNA qc barcode mapping */
int mk_DNA_qc_trgts(target * refs, opt_ref * ref_val, compact_target * ctrg, opt_BC * BC_val, association * assoc, FILE * fp_MUXtrgs, int trgt_ftype, target_params * trg_prms, int clcd_ctrg_cnt, TDSPLY_fasta * wt)
{
    extern char RLA29synch_3p11[34];
    
    int brcd_len = 0;  //barcode length
    int ctrg_cnt = 0;  //number of barcode targets;
    
    //test whether input targets file is in fasta format, which is required. if so, parse targets
    if (trgt_ftype == FASTA_FILE) {
        parse_DNA_qc_trgts(fp_MUXtrgs, trgt_ftype, ctrg, BC_val, assoc, trg_prms);
    } else {
        printf("mk_DNA_qc_trgts: error - targets must be provided in fasta format. aborting...\n");
        abort();
    }
    
    ctrg_cnt = trg_prms->t_cnt; //set ctrg_cnt, which tracks total ctrgs below
    
    printf("barcode count = %d\ncalc'd trg count = %d\n", trg_prms->t_cnt, clcd_ctrg_cnt);
    
    uint64_t i = 0;        //general purpose index
    uint64_t mutCode = 0;  //variable for storing mutation code. each barcode id is an unsigned 64 bit integer.
                           //bits 0-15 are used to store a mutation code. bits >= 16 are used to store the id
                           //of the input barcode from which the mutated barcode is derived
    
    int tot_SUB_trgts = 0; //number of substitution targets generated
    int tot_INS_trgts = 0; //number of insertion targets generated
    int tot_DEL_trgts = 0; //number of deletion targets generated
    
    for (i = 0; i < trg_prms->t_cnt; i++) { //for every input barcode
        
        mutCode = 1; //set mutCode to 1 (input barcodes have mutCode 0)
        
        //generate sub and indel targets for the current barcode
        tot_SUB_trgts += mk_SUB_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], trg_prms->BClen, &mutCode, 1);
        tot_INS_trgts += mk_INS_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], trg_prms->BClen, &mutCode, 1);
        tot_DEL_trgts += mk_DEL_trgts(ctrg, BC_val, &ctrg_cnt, &ctrg[i], trg_prms->BClen, &mutCode, toupper(RLA29synch_3p11[strlen(RLA29synch_3p11)-1]), 1);
    }
    
    //print the number of targets that were generated
    printf("total targets = %d\n", trg_prms->t_cnt + tot_SUB_trgts + tot_INS_trgts + tot_DEL_trgts);
    
    return ctrg_cnt;
}

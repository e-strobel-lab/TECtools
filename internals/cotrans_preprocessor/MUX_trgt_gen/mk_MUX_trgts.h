//
//  mk_MUX_trgts.h
//  
//
//  Created by Eric Strobel on 6/25/25.
//

#ifndef mk_MUX_trgts_h
#define mk_MUX_trgts_h

#include <stdio.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/constant_seqs.h"
#include "../../variant_maker/vmt_suffix.h"
#include "../../variant_maker/read_bcFile.h"

#include "../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/parse_fa_trgts.h"

#include "./mk_barcoded_target_fastas.h"

#define MIN_SUB_CODE 1    //minimum mutcode value for substitutions
#define MAX_SUB_CODE 48   //maximum mutcode value for substitutions
#define MIN_INS_CODE 49   //minimum mutcode value for insertions
#define MAX_INS_CODE 112  //maximum mutcode value for insertions
#define MIN_DEL_CODE 113  //minimum mutcode value for deletions
#define MAX_DEL_CODE 128  //maximum mutcode value for deletions

/* mk_MUX_trgts: manages TECprobe-MUX target generation */
int mk_MUX_trgts(TPROBE_names * nm, target * refs, opt_ref * ref_val, compact_target * ctrg, opt_BC * BC_val, FILE * fp_MUXtrgs, int trgt_ftype, target_params * trg_prms, int clcd_ctrg_cnt, TDSPLY_fasta * wt);

/* mk_SUB_trgts: generates single substitution targets for input barcode */
int mk_SUB_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode);

/* mk_INS_trgts: generate single insertion targets for input barcode */
int mk_INS_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode);

/* mk_DEL_trgts: generate single deletion targets for input barcode */
int mk_DEL_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode, char upstrm_nt);

/* get_target_type: determine target type using mutcode */
char * get_target_type(uint64_t mutcode, int * type_val);

#endif /* mk_MUX_trgts_h */

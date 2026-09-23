//
//  mk_DNA_qc_trgts.h
//  
//
//  Created by Eric Strobel on 9/9/26.
//

#ifndef mk_DNA_qc_trgts_h
#define mk_DNA_qc_trgts_h

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

/* mk_DNA_qc_trgts: generate targets for DNA qc barcode mapping */
int mk_DNA_qc_trgts(target * refs, opt_ref * ref_val, compact_target * ctrg, opt_BC * BC_val, association * assoc, FILE * fp_MUXtrgs, int trgt_ftype, target_params * trg_prms, int clcd_ctrg_cnt, TDSPLY_fasta * wt);

#endif /* mk_DNA_qc_trgts_h */

//
//  parse_fa_trgts.h
//  
//
//  Created by Eric Strobel on 10/7/25.
//

#ifndef parse_fa_trgts_h
#define parse_fa_trgts_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../cotrans_preprocessor/cotrans_preprocessor_defs.h"
#include "../../../cotrans_preprocessor/cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../../variant_maker/make_barcodes.h"
#include "../../../variant_maker/constant_seqs.h"

#include "./set_barcoded_compact_target.h"

/* parse_fa_trgts: stores input barcodes as targets */
void parse_fa_trgts(FILE * ifp, int trgt_ftype, void * trgts, void * trg_val, target_params * trg_prms);


#endif /* parse_fa_trgts_h */

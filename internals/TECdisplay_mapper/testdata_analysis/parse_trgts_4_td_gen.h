//
//  parse_trgts_4_td_gen.h
//  
//
//  Created by Eric Strobel on 7/9/26.
//

#ifndef parse_trgts_4_td_gen_h
#define parse_trgts_4_td_gen_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"

#include "../TECdisplay_mapper_defs.h"
#include "../TECdisplay_mapper_structs.h"

#include "../map_reads/map_expected/parse_vmt_trgts.h"

/* parse_trgts_4_td_gen: parse targets input file and store data in targets structures to be used for testdata generation. this function is used to circumvent the need to write a separate testdata generation function that is compatible with target sequences being stored in compact_target structures when using compact_target structures for read mapping */
void parse_trgts_4_td_gen(TDSPLY_names * nm, int trgt_ftype, target ** td_refs, opt_ref ** td_ref_val, target ** td_trgts, opt_mx_trg ** td_trg_val, target_params * td_trg_prms, TDSPLY_fasta * td_wt);

#endif /* parse_trgts_4_td_gen_h */

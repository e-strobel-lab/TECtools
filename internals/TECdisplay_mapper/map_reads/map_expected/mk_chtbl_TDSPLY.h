//
//  mk_chtbl_TDSPLY.h
//  
//
//  Created by Eric Strobel on 7/13/26.
//

#ifndef mk_chtbl_TDSPLY_h
#define mk_chtbl_TDSPLY_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

#include "../../../utils/debug.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "./get_key.h"

/* mk_chtbl_TDSPLY: construct hash table from target structure */
int mk_chtbl_TDSPLY(compact_h_node **chtbl, compact_h_node_bank *chn_bank, compact_target *ctrg, target *refs, target_params *trg_prms);

#endif /* mk_chtbl_TDSPLY_h */

//
//  mk_htbl_TDSPLY.h
//  
//
//  Created by Eric Strobel on 7/13/26.
//

#ifndef mk_htbl_TDSPLY_h
#define mk_htbl_TDSPLY_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

#include "../../../utils/debug.h"

#include "../../../seq_utils/seq2bin_hash.h"

#include "./get_key.h"

/* mk_htbl_TDSPLY: construct hash table from target structure */
int mk_htbl_TDSPLY(h_node **htbl, h_node_bank *bank, target *trgts, target *refs, target_params *trg_prms);



#endif /* mk_htbl_TDSPLY_h */

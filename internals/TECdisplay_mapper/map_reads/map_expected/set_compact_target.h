//
//  set_compact_target.h
//  
//
//  Created by Eric Strobel on 7/6/26.
//

#ifndef set_compact_target_h
#define set_compact_target_h

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

/* set_trgt: set target values in compact target struct */
void set_compact_target(compact_target * ctrg, opt_mx_trg * trg_val, target * crnt_ref, char * trgt_id, char * trgt_sq);

/* set_utility_target: set standard target as utility member of compact target struct */
void set_utility_target(compact_target * ctrg, target * utl, int cnt);

#endif /* set_compact_target_h */

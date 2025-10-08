//
//  set_trgt.h
//  
//
//  Created by Eric Strobel on 9/30/25.
//

#ifndef set_trgt_h
#define set_trgt_h

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
#include "./parse_vmt_trgts.h"

/* set_trgt: set target values in target struct */
void set_trgt(target * trgts, opt_mx_trg * trg_val, target * crnt_ref, char * trgt_id, char * trgt_sq);

#endif /* set_trgt_h */

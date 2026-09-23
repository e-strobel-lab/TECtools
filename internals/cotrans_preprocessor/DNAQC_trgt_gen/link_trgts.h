//
//  link_trgts.h
//  
//
//  Created by Eric Strobel on 9/10/26.
//

#ifndef link_trgts_h
#define link_trgts_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#define LINK_NTVS 0
#define LINK_MUTS 1

/* link_trgts: establish a link between two compact_targets using the 'lnk' member of the opt_BC structure within the compact targets. in LINK_NTVS mode, the lnk member of each opt_BC struct is pointed to the other compact target. in LINK_MUTS mode, a mutant target is provided alongside its parent target and the 'lnk' value of the mutant target is set to that of the parent target. */
void link_trgts(compact_target * trgt1, compact_target * trgt2, int mode);

#endif /* link_trgts_h */

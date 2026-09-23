//
//  link_trgts.c
//  
//
//  Created by Eric Strobel on 9/10/26.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "link_trgts.h"

/* link_trgts: establish a link between two compact_targets using the 'lnk' member of the opt_BC structure within the compact targets. in LINK_NTVS mode, the lnk member of each opt_BC struct is pointed to the other compact target. in LINK_MUTS mode, a mutant target is provided alongside its parent target and the 'lnk' value of the mutant target is set to that of the parent target. */
void link_trgts(compact_target * trgt1, compact_target * trgt2, int mode)
{
    opt_BC * BC_val1 = (opt_BC *)trgt1->opt; //optional barcode values for target 1
    opt_BC * BC_val2 = (opt_BC *)trgt2->opt; //optional barcode values for target 2
    
    if (mode == LINK_NTVS) {  //in LINK_NTVS mode
        BC_val1->lnk = trgt2; //point target 1's lnk to target 2
        BC_val2->lnk = trgt1; //point target 2's lnk to target 1
        
    } else if (mode == LINK_MUTS) {  //in LINK_MUTS mode
        BC_val2->lnk = BC_val1->lnk; //set target's lnk to the value of target 1's lnk
        
    } else { //unrecognized mode. throw error and abort
        printf("link_targets: error - unrecognized mode. aborting...\n");
        abort();
    }
    
    return;
}

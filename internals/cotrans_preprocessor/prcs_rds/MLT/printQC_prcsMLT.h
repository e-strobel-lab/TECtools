//
//  printQC_prcsMLT.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef printQC_prcsMLT_h
#define printQC_prcsMLT_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/gen_utils.h"
#include "testdata3pEnd_analysis.h"


/* print_metrics: prints a report of channel and 3' end processing metrics */
void print_metrics(names * nm, metrics  * met, fastp_params fastp_prms);

/* print_len_dist: print the observed 3' end distribution */
void print_len_dist(metrics *met, target3p_params trg_prms);

#endif /* printQC_prcsMLT_h */

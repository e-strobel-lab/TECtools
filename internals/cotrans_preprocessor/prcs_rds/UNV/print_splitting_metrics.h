//
//  print_splitting_metrics.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef print_splitting_metrics_h
#define print_splitting_metrics_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/gen_utils.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "../MLT/testdata3pEnd_analysis.h"
#include "../MUX/testdataMUX_analysis.h"

/* print_metrics: prints a report of channel and 3' end processing metrics */
void print_splitting_metrics(TPROBE_names * nm, mapping_metrics  * met, fastp_params fastp_prms);

/* print_len_dist: print the observed 3' end distribution */
void print_len_dist(mapping_metrics *met, target3p_params trg_prms);

#endif /* print_splitting_metrics_h */

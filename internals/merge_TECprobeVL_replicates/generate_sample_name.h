//
//  generate_sample_name.h
//  
//
//  Created by Eric Strobel on 1/26/24.
//

#ifndef generate_sample_name_h
#define generate_sample_name_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../utils/io_management.h"

#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"
#include "../cotrans_preprocessor/run_script_gen/MLT/mk_MLT_run_nm.h"

#include "parse_TECprobe_sample_name.h"

#include "../mkmtrx/cotrans_mtrx.h"
#include "../mkmtrx/mkmtrx_defs.h"

#include "./merge_TECprobeVL_replicates_defs.h"
#include "./merge_TECprobeVL_replicates_structs.h"

void generate_sample_name (sample_names * sn);
void merge_sample_names(sample_names * sn);

#endif /* generate_sample_name_h */

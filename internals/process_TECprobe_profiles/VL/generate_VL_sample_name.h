//
//  generate_VL_sample_name.h
//  
//
//  Created by Eric Strobel on 1/26/24.
//

#ifndef generate_VL_sample_name_h
#define generate_VL_sample_name_h

#include <stdio.h>

#include "../../global/global_defs.h"
#include "../../utils/io_management.h"

#include "../../cotrans_preprocessor/run_script_gen/UNV/config_struct.h"
#include "../../cotrans_preprocessor/run_script_gen/UNV/mk_run_nm.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"

#include "./parse_VL_sample_name.h"

/* generate_VL_sample_name: manages sample name parsing and generation */
void generate_VL_sample_name (sample_names * sn);


/* merge_sample_names: confirm that sample name attributes match and generate merged sample name */
void merge_sample_names(sample_names * sn);

#endif /* generate_VL_sample_name_h */

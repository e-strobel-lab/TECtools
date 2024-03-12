//
//  merge_VL_profiles.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef merge_VL_profiles_h
#define merge_VL_profiles_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"
#include "../global/store_SM2_profile.h" //TODO: temporary to allow compiling before revision of the function
#include "../global/calculate_normalization_factor.h"

/* merge_profiles: combine data from input reactivity profiles to generate merged reactivity profile files */
int merge_VL_profiles(SM2_analysis_directory * an_dir, int dir_count, SM2_analysis_directory * mrg, int min_depth, double max_bkg);



#endif /* merge_VL_profiles_h */

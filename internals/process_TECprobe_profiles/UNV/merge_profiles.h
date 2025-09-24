//
//  merge_profiles.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef merge_profiles_h
#define merge_profiles_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../../cotrans_preprocessor/cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor/run_script_gen/UNV/config_struct.h"

#include "../../seq_utils/isRNAbase.h"

#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"
#include "../UNV/calculate_normalization_factor.h"

/* merge_profiles: combine data from input reactivity profiles to generate merged reactivity profile files */
int merge_profiles(SM2_analysis_directory * an_dir, int dir_count, SM2_analysis_directory * mrg, int min_depth, double max_bkg);



#endif /* merge_profiles_h */

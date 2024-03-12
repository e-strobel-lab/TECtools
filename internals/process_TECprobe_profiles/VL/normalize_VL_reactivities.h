//
//  normalize_VL_reactivities.h
//  
//
//  Created by Eric Strobel on 2/16/24.
//

#ifndef normalize_VL_reactivities_h
#define normalize_VL_reactivities_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"
#include "../../utils/io_management.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../global/store_SM2_profile.h"
#include "../global/calculate_normalization_factor.h"

#include "./process_TECprobeVL_profiles_defs.h"
#include "./process_TECprobeVL_profiles_structs.h"


void normalize_VL_reactivities(SM2_analysis_directory * an_dir, int min_depth, double max_bkg, int verify);

#endif /* normalize_VL_reactivities_h */

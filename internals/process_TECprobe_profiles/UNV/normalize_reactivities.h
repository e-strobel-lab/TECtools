//
//  normalize_reactivities.h
//  
//
//  Created by Eric Strobel on 2/16/24.
//

#ifndef normalize_reactivities_h
#define normalize_reactivities_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"
#include "../../utils/io_management.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../UNV/store_SM2_profile.h"
#include "../UNV/calculate_normalization_factor.h"

#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"

/* normalize_reactivities: generate normalization factor using whole data set and normalize reactivity values */
void normalize_reactivities(SM2_analysis_directory * an_dir, int min_depth, double max_bkg, int norm_all, int verify_norm);

#endif /* normalize_reactivities_h */

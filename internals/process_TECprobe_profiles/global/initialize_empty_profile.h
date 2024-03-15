//
//  initialize_empty_profile.h
//  
//
//  Created by Eric Strobel on 3/15/24.
//

#ifndef initialize_empty_profile_h
#define initialize_empty_profile_h

#include <stdio.h>

#include "../../global/global_defs.h"
#include "../../utils/debug.h"
#include "../../utils/io_management.h"

#include "./store_SM2_profile.h"

/* initialize_empty_profile: initialize values for empty profile structure */
void initialize_empty_profile(SM2_profile * prf, int trg_nt_cnt, int trgt_start_ix);

#endif /* initialize_empty_profile_h */

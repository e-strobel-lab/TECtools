//
//  print_processed_profiles.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef print_processed_profiles_h
#define print_processed_profiles_h

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"

#include "../UNV/store_SM2_profile.h"

/* print_processed_profiles: generate output directories and files containing processed data */
void print_processed_profiles(SM2_analysis_directory * an_dir, char * out_dir, sample_names * sn);

/* print_profile: print shapemapper2 profile to file */
void print_profile(SM2_profile * prf, FILE * ofp);

#endif /* print_processed_profiles_h */

//
//  print_merged_profiles.h
//  
//
//  Created by Eric Strobel on 2/21/24.
//

#ifndef print_merged_profiles_h
#define print_merged_profiles_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"

#include "../UNV/store_SM2_profile.h"

/* print_merged_profiles: generated merged shapemapper2 output files */
void print_merged_profiles(SM2_analysis_directory * mrg, output_files * outfiles);

#endif /* print_merged_profiles_h */

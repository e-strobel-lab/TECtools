//
//  print_legacy_compiled_table.h
//  
//
//  Created by Eric Strobel on 2/24/24.
//

#ifndef print_legacy_compiled_table_h
#define print_legacy_compiled_table_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"

#include "../../mkmtrx/cotrans_mtrx.h"
#include "../../mkmtrx/mkmtrx_defs.h"

#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"

#include "../UNV/store_SM2_profile.h"

/* print_legacy_compiled_table: print aggregate TECprobe-VL profiles in the format used by Courtney's visualization tools */
void print_legacy_compiled_table(SM2_analysis_directory * an_dir, char * out_dir, sample_names * sn);

#endif /* print_legacy_compiled_table_h */

//
//  merge_reactivity_profiles.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef merge_reactivity_profiles_h
#define merge_reactivity_profiles_h

#include <stdio.h>
#include <stdlib.h>

#include "../global/global_defs.h"
#include "../mkmtrx/mkmtrx_defs.h"

#include "../cotrans_preprocessor/run_script_gen/MLT/config_MLT_struct.h"

#include "./merge_TECprobeVL_replicates_defs.h"
#include "./merge_TECprobeVL_replicates_structs.h"

/* merge_profiles: combine data from input reactivity profiles to generate merged reactivity profile files */
int merge_profiles(SM2_analysis_directory * an_dir, int dir_count, output_files * outfiles);

/* get_line_local: get line from file, place into array, remove trailing newline, and return
 line length if successful. local version that allows files to end on non-newline characters */
int get_line_local(char *line, FILE *ifp);

/* validate_header: confirm that column headers match expectations */
void validate_header(char * hdr);

#endif /* merge_reactivity_profiles_h */

//
//  merge_values_files.h
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#ifndef merge_values_files_h
#define merge_values_files_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

#include "./parse_TECdisplay_out_line.h"

/* merge_values_files: merge multiple values files into a single file */
int merge_values_files(values_input * vals_ipt, int vals_cnt, char * merged_out_nm, int nonstandard);

#endif /* merge_values_files_h */

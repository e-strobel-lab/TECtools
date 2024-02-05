//
//  make_output_directories.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef make_output_directories_h
#define make_output_directories_h

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "../global/global_defs.h"
#include "../mkmtrx/mkmtrx_defs.h"
#include "./merge_TECprobeVL_replicates_defs.h"
#include "./merge_TECprobeVL_replicates_structs.h"

/* make_output_directories: generate output directories and files */
void make_output_directories(SM2_analysis_directory * an_dir, output_files * outfiles, sample_names * sn);

#endif /* make_output_directories_h */

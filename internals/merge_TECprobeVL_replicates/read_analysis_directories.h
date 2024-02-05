//
//  read_analysis_directories.h
//  
//
//  Created by Eric Strobel on 1/25/24.
//

#ifndef read_analysis_directories_h
#define read_analysis_directories_h

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

#include "../global/global_defs.h"
#include "../mkmtrx/mkmtrx_defs.h"
#include "./merge_TECprobeVL_replicates_defs.h"
#include "./merge_TECprobeVL_replicates_structs.h"

/* read_prnt_directory: read parent shapemapper 2 analysis directory and identify transcript length analysis sub-folders*/
int read_prnt_directory(SM2_analysis_directory * an_dir, int dir_num, sample_names * sn);

/* read_tl_directory: read transcript length analysis directory and identify shapemapper 2 output sub-folders */
int read_tl_directory(SM2_analysis_directory * an_dir, int dir_num, int crnt_tl, char * tl_dir_nm, sample_names * sn);

/* read_SM2out_directory: read shapemapper 2 output directory and open reactivity profile file */
int read_SM2out_directory(SM2_analysis_directory * an_dir, int crnt_tl, char * out_dir_nm);

#endif /* read_analysis_directories_h */

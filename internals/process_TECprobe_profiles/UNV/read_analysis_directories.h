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
#include <limits.h>
#include <dirent.h>

#include "../../global/global_defs.h"
#include "../../mkmtrx/mkmtrx_defs.h"
#include "../process_TECprobe_profiles_defs.h"
#include "../process_TECprobe_profiles_structs.h"

#include "parse_sample_name.h"

/* read_prnt_directory: read parent shapemapper 2 analysis directory and identify target analysis directories */
int read_prnt_directory(SM2_analysis_directory * an_dir, int dir_num, sample_names * sn);

/* read_target_directory: read target directory and identify shapemapper 2 output sub-folders */
int read_target_directory(SM2_analysis_directory * an_dir, int dir_num, int crnt_tl, char * trg_dir_nm, sample_names * sn);

/* read_SM2out_directory: read shapemapper 2 output directory and open reactivity profile file */
int read_SM2out_directory(SM2_analysis_directory * an_dir, int crnt_tl, char * out_dir_nm);

/* get_nid: get numerical id of variant */
int get_nid(char * str, char * suffix);

/* set_len_range: set minimum and maximum target lengths */
void set_len_range(SM2_analysis_directory * an_dir);

/* test_trg_analysis_dir_format: test that target analysis directory
 name conforms to expected format. return length of id string. */
int test_trg_analysis_dir_format(char * str, char * suffix);

#endif /* read_analysis_directories_h */

//
//  process_output_files.h
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#ifndef process_output_files_h
#define process_output_files_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/stat.h>

#include "../utils/io_management.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

void clean_up_output(int layr_cnt, constraint_metadata * cons_meta);

void aggregate_output(int vals_cnt, values_input * vals, int layr_cnt, constraint_metadata * cons_meta, char * out_prefix, char * col_id);

void read_output_hdrs(int vals_cnt, char * col_id, int crnt_layr, FILE ** TECDnav_out, int * cols2incld, FILE ** ofp);

void read_data_line(int vals_cnt, int crnt_layr, FILE ** TECDnav_out, int * cols2incld, int * EOF_rchd, int * EOF_tot, FILE ** ofp);

#endif /* process_output_files_h */

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

/* clean_up_output: clean up output directory after TECdisplay_navigator analysis is complete */
void clean_up_output(int layr_cnt, constraint_metadata * cons_meta);

/* aggregate_output: aggregate all fracBound (or user-specified) columns for a given layer into an output file */
void aggregate_output(int vals_cnt, values_input * vals, int layr_cnt, constraint_metadata * cons_meta, char * out_prefix, char * col_id);

/* read_output_hdrs: read the header line of TECdisplay_navigator output files
 and print the headers of included columns to the aggregated data output file*/
void read_output_hdrs(int vals_cnt, char * col_id, int crnt_layr, FILE ** TECDnav_out, int * cols2incld, FILE ** ofp);

/* read_data_line: read data lines from TECdisplay_navigator output files and
 print the included columns to the relevant aggregated data output file*/
void read_data_line(int vals_cnt, int crnt_layr, FILE ** TECDnav_out, int * cols2incld, int * EOF_rchd, int * EOF_tot, FILE ** ofp);

#endif /* process_output_files_h */

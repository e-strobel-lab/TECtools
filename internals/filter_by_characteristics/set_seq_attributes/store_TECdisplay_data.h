//
//  store_TECdisplay_data.h
//  
//
//  Created by Eric Strobel on 3/5/26.
//

#ifndef store_TECdisplay_data_h
#define store_TECdisplay_data_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"
#include "../../utils/gen_utils.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../compare_variants/parse_comparison_file.h"

/* store_TECdisplay_data: store input TECdisplay data in the sequence attributes structure */
void store_TECdisplay_data(sequence_attributes * sq_att, char ** tecd_data_hdr, char *** prsd_tecd_data_hdr, int seq_cnt, char * tecd_path, comparison_values * cmp);

/* parse_tecd_data_hdr: parse TECdisplay data header string and store individual fields */
int parse_tecd_data_hdr(char *** prsd_tecd_data_hdr, char * tecd_data_hdr);

/* parse_tecd_data_line: parse TECdisplay data line and store values as doubles */
void parse_tecd_data_line(double ** vals, char * p_data, int fields);

/* find_cmp_index: identify index of comparison column specified by input comparison file */
void find_cmp_index(char ** prsd_tecd_data_hdr, int fields, comparison_values * cmp);

#endif /* store_TECdisplay_data_h */

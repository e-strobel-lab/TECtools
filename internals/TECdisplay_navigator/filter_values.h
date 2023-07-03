//
//  filter_values.h
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#ifndef filter_values_h
#define filter_values_h

#include <stdio.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"
#include "../seq_utils/ispair.h"
#include "../seq_utils/basemap.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"
#include "./read_vbase.h"
#include "./search_4_vbase_match.h"

/* filter_values: read input values file and send variants that match each set of
 constraints to the respective output file */
int filter_values(FILE * ipt, constraints * cons, int cons_cnt, basemap * bmap, char * out_dir_nm);

/* validate_data_line: validate data line from merged values file */
void validate_data_line(char * p_id, char * p_vals, int tot_vals);

/* test_bases: test if variant matches any input constraints */
int test_bases(char * p_id, base_params * vbases, int vbase_cnt, constraints * cons, int cons_cnt, basemap * bmap, int * match);

/* get_sample_info: parse column headers to obtain sample names and validate formatting */
int get_sample_info(char * line, sample_values * smpl_vals, int * tot_vals);

#endif /* filter_values_h */

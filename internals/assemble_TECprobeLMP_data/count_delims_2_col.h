//
//  count_delims_2_col.h
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#ifndef count_delims_2_col_h
#define count_delims_2_col_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "./assemble_TECprobeLMP_data_defs.h"
#include "./assemble_TECprobeLMP_data_structs.h"

int count_delims_2_col(FILE * ipt,  mode_parameters * mode_params, int nrchd_len, int * delims2col);

#endif /* count_delims_2_col_h */

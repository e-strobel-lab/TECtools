//
//  get_value.h
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#ifndef get_value_h
#define get_value_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

int get_value(FILE * ipt,  mode_parameters * mode_params, int delims2col, double * val);

#endif /* get_value_h */

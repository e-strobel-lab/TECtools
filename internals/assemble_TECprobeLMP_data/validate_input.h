//
//  validate_input.h
//  
//
//  Created by Eric Strobel on 2/8/24.
//

#ifndef validate_input_h
#define validate_input_h

#include <stdio.h>
#include <stdlib.h>

#include "../global/global_defs.h"
#include "../mkmtrx/cotrans_mtrx.h"
#include "../utils/io_management.h"

#include "./assemble_TECprobeLMP_data_defs.h"
#include "./assemble_TECprobeLMP_data_structs.h"

void validate_input (input_data * ipt);

#endif /* validate_input_h */

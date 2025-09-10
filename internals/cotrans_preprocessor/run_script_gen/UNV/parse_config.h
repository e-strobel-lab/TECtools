//
//  parse_config.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef parse_config_h
#define parse_config_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/gen_utils.h"
#include "../../../utils/io_management.h"
#include "config_struct.h"

/* parse_config: set config_MLT structure variables to config file settings */
int parse_config(FILE * ifp, tprobe_configuration * config, int * mode);

#endif /* parse_config_h */

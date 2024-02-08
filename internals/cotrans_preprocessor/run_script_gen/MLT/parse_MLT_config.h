//
//  parse_MLT_config.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef parse_MLT_config_h
#define parse_MLT_config_h

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
#include "config_MLT_struct.h"

/* parse_MLT_config: set config_MLT structure variables to config file settings */
int parse_MLT_config(FILE * ifp, configuration_MLT * config_MLT);

#endif /* parse_MLT_config_h */

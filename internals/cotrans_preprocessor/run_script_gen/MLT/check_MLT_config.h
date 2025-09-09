//
//  check_MLT_config.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef check_MLT_config_h
#define check_MLT_config_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "config_MLT_struct.h"

//check_MLT_config: check that config_MLT variables have been set properly
int check_MLT_config(configuration_MLT * config_MLT, int mode);

#endif /* check_MLT_config_h */

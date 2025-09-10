//
//  check_config.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef check_config_h
#define check_config_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "config_struct.h"

//check_config: check that config_MLT variables have been set properly
int check_config(tprobe_configuration * config, int mode);

#endif /* check_config_h */

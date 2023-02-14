//
//  mk_MLT_config.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef mk_MLT_config_h
#define mk_MLT_config_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

/* mk_MLT_config:
 generate configuration file for use make_shapemapper2_run_script executable
 this config file is used to generate a shell script that run commands for
 shapemapper2 analysis */
void mk_MLT_config(names * nm, target3p_params trg_prms);

#endif /* mk_MLT_config_h */

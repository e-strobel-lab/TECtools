//
//  mk_config.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef mk_config_h
#define mk_config_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "../../run_script_gen/UNV/mk_run_script.h"

/* mk_config:
 generate configuration file for use make_shapemapper2_run_script executable
 this config file is used to generate a shell script that run commands for
 shapemapper2 analysis */
void mk_config(TPROBE_names * nm, target3p_params * trg_prms, int mode);

#endif /* mk_config_h */

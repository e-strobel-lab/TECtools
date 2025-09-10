//
//  mk_run_script.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef mk_run_script_h
#define mk_run_script_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/debug.h"
#include "config_struct.h"
#include "parse_config.h"
#include "check_config.h"
#include "../MLT/mk_MLT_run_nm.h"
#include "../MLT/print_MLT_SM2_script.h"


/* mk_run_script: manages shapemapper2 run script generation for multilength cotranscriptional RNA
 structure probing experiments */
int mk_run_script(FILE * fp_config);

/* check initialization: check initialization of config struct */
void check_cnfgMLT_init(tprobe_configuration * config);

#endif /* mk_run_script_h */

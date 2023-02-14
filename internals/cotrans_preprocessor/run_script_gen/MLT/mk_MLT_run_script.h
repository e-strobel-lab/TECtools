//
//  mk_MLT_run_script.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef mk_MLT_run_script_h
#define mk_MLT_run_script_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/debug.h"
#include "config_MLT_struct.h"
#include "parse_MLT_config.h"
#include "check_MLT_config.h"
#include "mk_MLT_run_nm.h"
#include "print_MLT_SM2_script.h"


/* mk_MLT_run_script: manages shapemapper2 run script generation for multilength cotranscriptional RNA
 structure probing experiments */
int mk_MLT_run_script(FILE * fp_config_MLT);

/* check initialization: check initialization of config struct */
void check_cnfgMLT_init(configuration_MLT * config_MLT);

#endif /* mk_MLT_run_script_h */

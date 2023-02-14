//
//  mk_MLT_run_nm.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef mk_MLT_run_nm_h
#define mk_MLT_run_nm_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "config_MLT_struct.h"

//mk_MLT_run_nm: generate sample name from config file settings
int mk_MLT_run_nm(char * sample_name, configuration_MLT * config_MLT);

#endif /* mk_MLT_run_nm_h */

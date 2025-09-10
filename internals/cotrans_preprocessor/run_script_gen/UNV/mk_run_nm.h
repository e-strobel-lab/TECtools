//
//  mk_run_nm.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef mk_run_nm_h
#define mk_run_nm_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../UNV/config_struct.h"

//mk_run_nm: generate sample name from config file settings
int mk_run_nm(char * sample_name, tprobe_configuration * config);

#endif /* mk_run_nm_h */

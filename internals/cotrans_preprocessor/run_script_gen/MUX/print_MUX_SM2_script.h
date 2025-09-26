//
//  print_MUX_SM2_script.h
//  
//
//  Created by Eric Strobel on 9/24/25.
//

#ifndef print_MUX_SM2_script_h
#define print_MUX_SM2_script_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../UNV/config_struct.h"

/* print_MUX_SM2_script: print shell scripts that run shapemapper2 for all barcoded targets */
int print_MUX_SM2_script(char * sample_name, int brcd_cnt, char ** brcd_id, tprobe_configuration * config);

#endif /* print_MUX_SM2_script_h */

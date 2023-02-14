//
//  print_MLT_SM2_script.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef print_MLT_SM2_script_h
#define print_MLT_SM2_script_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "config_MLT_struct.h"

/*print_MLT_SM2_script: print shell script that runs shapemapper2 for all 3' end lengths */
int print_MLT_SM2_script(char * sample_name, configuration_MLT * config_MLT);

#endif /* print_MLT_SM2_script_h */

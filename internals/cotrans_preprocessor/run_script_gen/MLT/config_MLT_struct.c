//
//  config_MLT_struct.c
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "config_MLT_struct.h"

/* set_cfg_string: allocate memory for configuration_MLT structure string */
void set_cfg_string(char ** cfg_str, char * val, int buffer) {
    
    if (((*cfg_str) = calloc(strlen(val)+buffer+1, sizeof(**cfg_str))) == NULL) {
        printf("alloc_cfg_str_mem: memory allocation for configuration_MLT string failed. aborting...\n");
        abort();
    }
    
    strcpy(*cfg_str, val);
}

/* set_TF_value: set true or false value as 1 and 0, respectively */
void set_TF_value(char * tf, char * setting_name, int * config_val)
{
    if (!strcmp(tf, "FALSE")) {
        *config_val = 0;
    } else if (!strcmp(tf, "TRUE")) {
        *config_val = 1;
    } else {
        printf("parse_MLT_config: error - unexpected value for %s setting. set to FALSE, or TRUE. aborting...\n", setting_name);
        abort();
    }
}

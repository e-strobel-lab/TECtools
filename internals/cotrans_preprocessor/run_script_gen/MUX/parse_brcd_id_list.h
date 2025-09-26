//
//  parse_brcd_id_list.h
//  
//
//  Created by Eric Strobel on 9/24/25.
//

#ifndef parse_brcd_id_list_h
#define parse_brcd_id_list_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"

#include "../UNV/config_struct.h"

/* parse_brcd_id_list: parse barcode ids file and store barcode ids */
int parse_brcd_id_list(tprobe_configuration * config, char *** brcd_id);

#endif /* parse_brcd_id_list_h */

//
//  get_constraint_metadata.h
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#ifndef get_constraint_metadata_h
#define get_constraint_metadata_h

#include <stdio.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"

#include "../seq_utils/basemap.h"

#include "../TECdisplay_navigator/parse_reference.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

void get_constraint_metadata(char * ipt_fn, constraint_metadata * cons_meta, char type);

#endif /* get_constraint_metadata_h */

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

/* get_constraint_metadata: read constraint files to get information about individual constraints within
 the file. This process also validates each constraint before it is used by TECdisplay_navigator */
void get_constraint_metadata(char * ipt_fn, constraint_metadata * cons_meta, char type, int layr_cnt);

/* valid8_constraint_compatiblity: confirm all constraints files are compatible */
void valid8_constraint_compatiblity(int layr_cnt, constraint_metadata * cons_meta);

/* print_constraint_metadata: print metadata for each input constraint file to output file */
void print_constraint_metadata(int layr_cnt, constraint_metadata * cons_meta);

#endif /* get_constraint_metadata_h */

//
//  Hfilter.h
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#ifndef Hfilter_h
#define Hfilter_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "../utils/io_management.h"

#include "../seq_utils/basemap.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

void Hfilter(constraint_metadata * cons_meta, char * prev_out_prfx, int cl, int layr_cnt, char * prev_TECDnav_out, char * layr_list[MAX_LAYERS]);

void store_out_filenames(int cl, char * prefix, constraint_metadata * cons_meta);

#endif /* Hfilter_h */

//
//  calc_output_files.h
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#ifndef calc_output_files_h
#define calc_output_files_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

/* calc_output_files: calculate the expected number of output files and
 confirm that TECdisplay_navigator analysis should proceed */
int calc_output_files(int layr_cnt, constraint_metadata * cons_meta);

#endif /* calc_output_files_h */

//
//  TECdisplay_Hnav_global_vars.c
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#include <stdio.h>

#include "../global/global_defs.h"

#include "../TECdisplay_navigator/TECdisplay_navigator_defs.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

#include "TECdisplay_Hnav_global_vars.h"

//global variables
char TDHN_TECDnav_path[MAX_LINE] = {0};            //TECdisplay_navigator path
char TDHN_merged_out_nm[15] = {"merged_out.txt"};  //merged_output filename

output_file_names out_fns[MAX_LAYERS] = {0};

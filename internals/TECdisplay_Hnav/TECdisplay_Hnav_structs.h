//
//  TECdisplay_Hnav_structs.h
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#ifndef TECdisplay_Hnav_structs_h
#define TECdisplay_Hnav_structs_h

#include <stdio.h>

#include "../global/global_defs.h"

#include "../seq_utils/basemap.h"

#include "../TECdisplay_navigator/TECdisplay_navigator_defs.h"
#include "../TECdisplay_navigator/TECdisplay_navigator_structs.h"

#include "./TECdisplay_Hnav_defs.h"



typedef struct constraint_metadata {
    FILE * fp;
    char fn[MAX_NAME];
    char sn[MAX_NAME];
    char cn[MAX_CONSTRAINTS][MAX_NAME];
    char typ;
    wt_source wt;
    basemap bmap;
    char * cnstnt_indels;
    int c_cnt;
} constraint_metadata;

typedef struct output_file_names {
    char ** fn;
    int f_cnt;
    int nxt;
} output_file_names;

#endif /* TECdisplay_Hnav_structs_h */

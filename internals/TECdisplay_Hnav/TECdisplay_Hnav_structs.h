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


/* constraint_metadata: storage for constraint file information*/
typedef struct constraint_metadata {
    FILE * fp;                               //pointer to constraint file
    char fn[MAX_NAME];                       //constraint filename
    char sn[MAX_NAME];                       //constraint sample name (no file suffix)
    char * cn[MAX_CONSTRAINTS];              //constraint names
    char typ;                                //constraint type: c = included constraint, x = excluded constraint
    wt_source wt;                            //wt source sequence
    basemap bmap;                            //basemap of variant library
    char * cnstnt_indels;                    //constant indels in the variant library
    struct constraints con[MAX_CONSTRAINTS]; //storage for variant constraints
    int c_cnt;                               //number of contraints in constraint file
} constraint_metadata;

/* output_file_names: storage for output file names*/
typedef struct output_file_names {
    char ** fn;      //pointers for allocating file name storage
    int f_cnt;       //number of output file names
    int nxt;         //next output file name index
} output_file_names;

#endif /* TECdisplay_Hnav_structs_h */

//
//  mapping_metrics.c
//  
//
//  Created by Eric Strobel on 10/14/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "mapping_metrics.h"

/* init_chnl_mtrcs_mem: initialize channel tracking metrics memory */
void init_chnl_mtrcs_mem(mapping_metrics * met, int channels)
{
    if ((met->chan_count = calloc(channels, sizeof(*(met->chan_count)))) == NULL) {
        printf("init_chnl_mtrcs_mem: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    if ((met->full_match = calloc(channels, sizeof(*(met->full_match)))) == NULL) {
        printf("init_chnl_mtrcs_mem: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    if ((met->part_match = calloc(channels, sizeof(*(met->part_match)))) == NULL) {
        printf("init_chnl_mtrcs_mem: error - memory allocation failed. aborting...\n");
        abort();
    }
    
    return;
}

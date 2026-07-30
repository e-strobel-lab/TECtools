//
//  store_TECdisplay_data.h
//  
//
//  Created by Eric Strobel on 3/5/26.
//

#ifndef store_TECdisplay_data_h
#define store_TECdisplay_data_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"
#include "../../utils/gen_utils.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

/* store_TECdisplay_data: store input TECdisplay data in the sequence attributes structure */
void store_TECdisplay_data(sequence_attributes * sq_att, char ** tecd_data_hdr, int seq_cnt, char * tecd_path);

#endif /* store_TECdisplay_data_h */

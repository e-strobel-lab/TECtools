//
//  init_seq_attributes.h
//  
//
//  Created by Eric Strobel on 11/21/25.
//

#ifndef init_seq_attributes_h
#define init_seq_attributes_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

/* init_seq_attributes: initialize sequence attributes structure */
void init_seq_attributes(sequence_attributes ** sq_att, descriptor * des, descriptor_bank * des_bnk, int seq_cnt, int des_cnt, int * typ_cnt);

#endif /* init_seq_attributes_h */

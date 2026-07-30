//
//  print_annotated_data.h
//  
//
//  Created by Eric Strobel on 3/5/26.
//

#ifndef print_annotated_data_h
#define print_annotated_data_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

#include "../set_seq_attributes/set_attributes.h"


/* print_annotated_data: print output file that includes all sequences, associated data, and structure predictions */
void print_annotated_data(sequence_attributes * sq_att, int seq_cnt, int des_cnt, char * tecd_data_hdr);

#endif /* print_annotated_data_h */

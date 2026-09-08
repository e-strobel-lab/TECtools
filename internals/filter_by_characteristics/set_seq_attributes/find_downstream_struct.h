//
//  find_downstream_struct.h
//  
//
//  Created by Eric Strobel on 8/25/26.
//

#ifndef find_downstream_struct_h
#define find_downstream_struct_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../seq_utils/mk_fasta.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

#include "./parse_descriptor_file.h"
#include "./get_msa_subseq.h"
#include "./run_RNAStructure.h"
#include "./parse_ct_file.h"
#include "./set_attributes.h"

//NOTE: need control on whether this is run - don't need to run extra predictions if assessing full sequence, prob not necessary with reactivity constraints either
//NOTE: currently only works for proximal deltaG predictions
/* find_downstream_struct: perform iterative structure predictions to identify structures downstream of an initial structure */
int find_downstream_struct(structProps * sp, descriptor * des, sequence_attributes * sq_att, char * path2RNAStructure, char * nm, char * path2ct, int ct_path_maxlen, int ptyp, int ext_limit);

/* mv_ct_to_LL_top: move a con_table struct to the top of the con_table linked list (not currently using this) */
void mv_ct_to_LL_top(con_table * new_top, con_table * old_top);

#endif /* find_downstream_struct_h */

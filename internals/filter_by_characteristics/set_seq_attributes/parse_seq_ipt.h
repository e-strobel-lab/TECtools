//
//  parse_seq_ipt.h
//  
//
//  Created by Eric Strobel on 11/20/25.
//

#ifndef parse_seq_ipt_h
#define parse_seq_ipt_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/isDNAbase.h"
#include "../../seq_utils/isRNAbase.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

/* parse_seq_ipt: parse input multiple sequence alignment file */
//note: aligment characters are preseved in the stored sequences and removed later as needed when generating descriptor sequences
int parse_seq_ipt(char * ipt_path, int msa_typ, sequence_attributes ** sq_att);

#endif /* parse_seq_ipt_h */

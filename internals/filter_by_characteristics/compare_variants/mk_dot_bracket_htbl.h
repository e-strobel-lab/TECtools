//
//  mk_dot_bracket_htbl.h
//  
//
//  Created by Eric Strobel on 6/30/26.
//

#ifndef mk_dot_bracket_htbl_h
#define mk_dot_bracket_htbl_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../set_seq_attributes/parse_descriptor_file.h"

/* opt_db: optional values struct for associating structProps structure with a compact_target structure */
typedef struct opt_db {
    structProps ** sp;
} opt_db;

/* mk_dot_bracket_htbl: makes compact target hash table */
/* hash table has linked list buckets for possible collisions */
void mk_dot_bracket_htbl(compact_h_node ** htbl, compact_h_node_bank * bank, compact_target * ctrg, int ctrg_cnt, structProps ** sp_list/*, target_params * trg_prms*/);

/* dot_bracket_to_bin_hash: convert DNA sequence to two bit notation */
int dot_bracket_to_bin_hash(char * hash_seq, binary_seq * bin_seq, int array_max);

/* init_opt_db: allocate opt_db structure memory in compact_target structures */
void init_opt_db(compact_target * ctrg, int ctrg_cnt);

/* map_structProps_2_htbl: associate structProps structures with the relevant targets */
void map_structProps_2_htbl(structProps ** sp_list, int sp_cnt, compact_h_node ** htbl);


#endif /* mk_dot_bracket_htbl_h */

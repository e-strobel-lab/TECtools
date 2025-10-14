//
//  map_barcoded_targets.h
//  
//
//  Created by Eric Strobel on 10/14/25.
//

#ifndef map_barcoded_targets_h
#define map_barcoded_targets_h

#include <stdio.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../../utils/io_management.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "../../MUX_trgt_gen/mk_MUX_trgts.h"
#include "../../MUX_trgt_gen/set_barcoded_compact_target.h"

/* mk_htbl_MUX: makes compact target hash table */
/* hash table has linked list buckets for possible collisions */
void mk_htbl_MUX(compact_h_node ** htbl_MUX, compact_h_node_bank * bank, compact_target * ctrg, int count, mapping_metrics * met);

/* hash_brcd_trgt: generates hash key for binary encoded sequence (up to 32 nt/64 bits */
uint64_t hash_brcd_trgt(binary_seq * bsq);

/* check_brcd_diff: checks whether redundant barcode targets are linked to different barcode variants*/
int check_brcd_diff(compact_target * old, compact_target * new);

/* map_brcd: map barcode to target using hash table */
compact_target * map_brcd(char * brcd_str, char * rc_brcd_str, compact_h_node **htbl_MUX, compact_target ** mpd_trg, mapping_metrics * met);

/* get_target_type: determine target type using mutcode */
char * get_target_type(uint64_t mutcode, int * type_val);

#endif /* map_barcoded_targets_h */

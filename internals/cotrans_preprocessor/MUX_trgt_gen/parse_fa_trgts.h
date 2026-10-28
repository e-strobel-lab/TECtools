//
//  parse_fa_trgts.h
//  
//
//  Created by Eric Strobel on 10/7/25.
//

#ifndef parse_fa_trgts_h
#define parse_fa_trgts_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"
#include "../../seq_utils/C3SCbc_tree.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"
#include "../../variant_maker/bc_ind.h"

#include "./parse_id_link.h"
#include "./set_barcoded_compact_target.h"

#define ID_TYPE_INIT 0
#define CONCISE_VAR_ID 1
#define VERBOSE_VAR_ID 2
#define NONSTD_VAR_ID  3

#define MAX_VAR_NUM 99999
#define MAX_NUM_STR 10

/* parse_fa_trgts: stores input barcodes as targets */
void parse_fa_trgts(FILE * ifp, int trgt_ftype, int lib_struct, void * trgts, void * trg_val, C3SC_barcode_node * scbt_root, var_id_link * id_link, target_params * trg_prms, int data_type);

/* detect_id_type: determine variant id type */
int detect_id_type(char * ipt_trgt_id);

/* sub_vrbs_id: use var_id_link structure to substitute concise id with verbose id */
void sub_vrbs_id(char ** id2use, char * ipt_trgt_id, var_id_link * id_link, int cnt);

/* append var_num: append variant number to id. for use with non-standard variant ids */
void append_var_num(char ** id2use, char * ipt_trgt_id);

#endif /* parse_fa_trgts_h */

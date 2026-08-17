//
//  cmpr_smlr_structs.h
//  
//
//  Created by Eric Strobel on 6/29/26.
//

#ifndef cmpr_smlr_structs_h
#define cmpr_smlr_structs_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/data_dist.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../set_seq_attributes/parse_descriptor_file.h"

#include "../set_seq_attributes/parse_ct_file.h"

#include "./parse_comparison_file.h"
#include "./mk_dot_bracket_htbl.h"

/* cmpr_smlr_structs: perform a systematic comparison of variants with closely-related structures to identify variants that exhibit large effects on function */
void cmpr_smlr_structs(sequence_attributes * sq_att, descriptor * des, int seq_cnt, int des_cnt, char * des_nm, comparison_values * cmp);

/* set_db_compact_target: generate compact_targets from structProps structures and associate the structProps struct with the compact_target */
void set_db_compact_target(compact_target * ctrg, opt_db * db_val, structProps ** sp_list, int ctrg_cnt, sequence_attributes * sq_att, int seq_cnt, int des_ix, int des_typ);

/* set_sp_z_scores: set z-scores for comparison values */
void set_sp_z_scores(compact_target * ctrg, int ctrg_cnt, comparison_values * cmp);

/* find_ntbl_var_prs: find variant pairs that exhibit closely related structures but substantial functional differences */
void find_ntbl_var_prs(compact_target * ctrg, int ctrg_cnt, comparison_values * cmp);

/* sp_cmpfunc: simple float comparison function used to sort structProps array. */
int sp_cmpfnc(const void * a, const void * b);

/* cnt_dif_bps: count differences between two related structures */
void cnt_dif_bps(mct_diffs * dif, min_con_table * mct1, min_con_table * mct2);

#endif /* cmpr_smlr_structs_h */

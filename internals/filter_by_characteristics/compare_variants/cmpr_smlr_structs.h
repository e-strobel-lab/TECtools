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
#include <math.h>
#include <float.h>
#include <unistd.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../utils/data_dist.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../set_seq_attributes/parse_descriptor_file.h"

#include "../set_seq_attributes/parse_ct_file.h"

#include "./mk_dot_bracket_htbl.h"
#include "./parse_comparison_file.h"

#define SIM 0 //similar structures mode code
#define DIF 1 //different structures mode code

#define MIN_Z_DIF_INIT  2.0 //minimum z-score difference, used when looking for functional differences
#define MAX_Z_DIF_INIT  0.5 //maximum z-score difference, used when looking for funcitonal similarities
#define MAX_DG_DIF_INIT 0.5 //maximum delta G difference


/* cmpr_smlr_structs: perform a systematic comparison of variants with closely-related structures to identify variants that exhibit large effects on function */
void cmpr_smlr_structs(sequence_attributes * sq_att, descriptor * des, int seq_cnt, int des_cnt, comparison_values * cmp);

/* set_db_compact_target: generate compact_targets from structProps structures and associate the structProps struct with the compact_target */
void set_db_compact_target(compact_target * ctrg, opt_db * db_val, structProps ** sp_list, int ctrg_cnt, sequence_attributes * sq_att, int seq_cnt, int des_ix, int des_typ);

/* ctrg_cmp_fnc: comparison function for sorting pointers to compact_target pointers by count */
int ctrg_cmpfnc(const void * a, const void * b);

/* print_structure_inventory: prints an inventory of identified structure classes */
void print_structure_inventory(compact_target ** ctrg, int wl, comparison_values * cmp, char * dir_nm);

/* set_sp_z_scores: set z-scores for comparison values */
void set_sp_z_scores(compact_target ** ctrg, int wl, comparison_values * cmp);

/* find_ntbl_var_prs: find variant pairs that exhibit closely related structures but substantial functional differences */
void find_ntbl_var_prs(compact_target ** ctrg, int wl, comparison_values * cmp, char * dir_nm);

/* sp_cmpfunc: simple float comparison function used to sort structProps array. */
int sp_cmpfnc(const void * a, const void * b);

/* cnt_dif_bps: count differences between two related structures */
//NOTE: cssq memory needs to be freed outside function
void cnt_dif_bps(mct_diffs * difs, min_con_table * mct1, min_con_table * mct2);

/* print_tl_out_hdr: print header line for two-line output file */
void print_tl_out_hdr(FILE * ofp, char * db);

/* print_tl_out_data: print data lines for two-line output file */
void print_tl_out_data(FILE * ofp, comparison_values * cmp, sequence_attributes * sq_att1, sequence_attributes * sq_att2, structProps * sp1, structProps * sp2, mct_diffs * difs);

#endif /* cmpr_smlr_structs_h */

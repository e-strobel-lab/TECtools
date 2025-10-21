//
//  map_standard_TDSPLY_reads.h
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#ifndef map_standard_TDSPLY_reads_h
#define map_standard_TDSPLY_reads_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../TECdisplay_mapper_defs.h"
#include "../TECdisplay_mapper_structs.h"

#include "../../utils/io_management.h"
#include "../../utils/debug.h"
#include "../../seq_utils/seq2bin_hash.h"
#include "./UNV/call_fastp_TDSPLY.h"
#include "./UNV/prcs_chnl_TDSPLY.h"
#include "./map_expected/get_key.h"
#include "../../seq_utils/revcomp.h"
#include "../../seq_utils/basemap.h"
#include "../../seq_utils/mapping_metrics.h"
#include "./map_expected/parse_vmt_trgts.h"
#include "./map_expected/mk_output_files.h"
#include "./map_expected/print_navigator_template.h"
#include "../testdata_analysis/mk_TDSPLY_test_data.h"
#include "../testdata_analysis/assess_TDSPLY_test_data.h"

/* prcs_standard_TDSPLY_reads: coordinates targets parsing, fastp processing, read mapping, and output file gen */
int prcs_standard_TDSPLY_reads(TDSPLY_names * nm, FILE * fp_trgs, int trgt_ftype, char * minQ, fastp_params fastp_prms, testdata_vars * testdata, int mode);

/* mk_htbl_TDSPLY: construct hash table from target structure */
int mk_htbl_TDSPLY(h_node **htbl, h_node_bank *bank, target *trgts, target *refs, target_params *trg_prms);

/* map_standard_TDSPLY_reads: map reads to user-supplied vmt targets */
void map_standard_TDSPLY_reads(FILE *ifp, h_node **htbl, target *refs, target *trgts, char * minQ, target_params * trg_prms, mapping_metrics * met, testdata_vars * testdata, int mode);

/* test_cbase_qscores: test whether all constant bases meet or exceed the minimum qscore */
int test_cbase_qscores(char * qscore5p, char * minQc, target *refs);

/* count_matched_targets: count the number of targets to which at least 1 read mapped */
int count_matched_targets(target * trgts, target_params * trg_prms);

/* crrct_testdata_nonsrc_mtch: decrement match counters when a mutant testdata read maps to a target
 from a different source sequence. this allows testdata analysis to be run correctly using targets
 that were generated from very closely related variant templates */
void crrct_testdata_nonsrc_mtch(target * trg, opt_mx_trg * trg_vals, mapping_metrics * met, int channel, int chnl_mtch_typ, testdata_vars * testdata);

#endif /* map_standard_TDSPLY_reads_h */

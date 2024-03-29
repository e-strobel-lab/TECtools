//
//  map_reads.h
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#ifndef map_reads_h
#define map_reads_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../TECdisplay_mapper_defs.h"
#include "../TECdisplay_mapper_structs.h"

#include "../../utils/io_management.h"
#include "../../utils/debug.h"
#include "../../seq_utils/seq2bin_hash.h"
#include "./UNV/call_fastp.h"
#include "../../seq_utils/revcomp.h"
#include "../../seq_utils/basemap.h"
#include "./map_expected/parse_mx_trgts.h"
#include "./map_expected/mk_output_files.h"
#include "../testdata_analysis/mk_test_data.h"
#include "../testdata_analysis/assess_test_data.h"

/* map_reads: coordinates targets parsing, fastp processing, read mapping, and output file generation */
int map_reads (names * nm, FILE * fp_trgs, char * minQ, fastp_params fastp_prms, testdata_vars * testdata, int mode);

/* mk_htbl_TDSPLY: construct hash table from target structure */
int mk_htbl_TDSPLY(h_node **htbl, h_node_bank *bank, target *trgts, target *refs, target_params *trg_prms);

/* map_expected_reads: map reads to user-supplied targets */
void map_expected_reads(FILE *ifp, h_node **htbl, target *refs, target *trgts, char * minQ, target_params * trg_prms, metrics * met, testdata_vars * testdata, int mode);

/* get_key: generate key string composed the nucleotides at variable
 base positions in the input read sequence */
int get_key(char * key, char * end5p, char * qscore5p, char * minQv, target *refs, int key_type);

/* test_cbase_qscores: test whether all constant bases meet or exceed the minimum qscore */
int test_cbase_qscores(char * qscore5p, char * minQc, target *refs);

/* count_matched_targets: count the number of targets to which at least 1 read mapped */
int count_matched_targets(target * trgts, target_params * trg_prms);

/* crrct_testdata_nonsrc_mtch: decrement match counters when a mutant testdata read maps to a target
 from a different source sequence. this allows testdata analysis to be run correctly using targets
 that were generated from very closely related variant templates */
void crrct_testdata_nonsrc_mtch(target * trg, opt_mx_trg * trg_vals, metrics * met, int channel, int chnl_mtch_typ, testdata_vars * testdata);



#endif /* map_reads_h */

//
//  map_barcoded_TDSPLY_reads.h
//  
//
//  Created by Eric Strobel on 10/14/25.
//

#ifndef map_barcoded_TDSPLY_reads_h
#define map_barcoded_TDSPLY_reads_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../cotrans_preprocessor/MUX_trgt_gen/mk_MUX_trgts.h"
#include "../../cotrans_preprocessor/prcs_rds/MUX/map_barcoded_targets.h"
#include "../../cotrans_preprocessor/prcs_rds/MUX/get_brcd_str.h"
#include "../../cotrans_preprocessor/prcs_rds/MUX/testdataMUX_analysis.h"

#include "../TECdisplay_mapper_defs.h"
#include "../TECdisplay_mapper_structs.h"

#include "../../utils/io_management.h"
#include "../../utils/debug.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "./UNV/call_fastp_TDSPLY.h"
#include "./UNV/prcs_chnl_TDSPLY.h"

#include "../../seq_utils/revcomp.h"
#include "../../seq_utils/basemap.h"
#include "../../seq_utils/mapping_metrics.h"

#include "./map_expected/mk_output_files.h"
#include "./map_expected/print_navigator_template.h"

#include "../testdata_analysis/mk_TDSPLY_test_data.h"
#include "../testdata_analysis/assess_TDSPLY_test_data.h"

/* prcs_barcoded_TDSPLY_reads: coordinates targets parsing, fastp processing, read mapping, and output file generation */
int prcs_barcoded_TDSPLY_reads(TDSPLY_names * nm, FILE * fp_trgs, int trgt_ftype, char * minQ, fastp_params fastp_prms, testdata_vars * testdata, int mode);

/* map_barcoded_TDSPLY_reads: map reads to user-supplied fasta barcoded targets */
void map_barcoded_TDSPLY_reads(FILE *ifp, compact_h_node **htbl, compact_target *ctrg, target_params * trg_prms, mapping_metrics * met, int mode);

/* count_matched_compact_targets: count the number of copmact targets to which at least 1 read mapped */
int count_matched_compact_targets(compact_target * ctrg, target_params * trg_prms, int ctrg_cnt);

#endif /* map_barcoded_TDSPLY_reads_h */

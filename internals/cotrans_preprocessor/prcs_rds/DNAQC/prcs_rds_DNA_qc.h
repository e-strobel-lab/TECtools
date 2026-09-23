//
//  prcs_rds_DNA_qc.h
//  
//
//  Created by Eric Strobel on 9/8/26.
//

#ifndef prcs_rds_DNA_qc_h
#define prcs_rds_DNA_qc_h

#include <stdio.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../../variant_maker/variant_maker_defs.h"
#include "../../../variant_maker/vmt_suffix.h"
#include "../../../variant_maker/make_barcodes.h"

#include "../../../TECdisplay_mapper/TECdisplay_mapper_defs.h"
#include "../../../TECdisplay_mapper/TECdisplay_mapper_structs.h"
#include "../../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"
#include "../../../TECdisplay_mapper/map_reads/map_expected/print_navigator_template.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "../UNV/bypass_fastp.h"
#include "../MUX/map_barcoded_targets.h"

#include "../../DNAQC_trgt_gen/mk_DNA_qc_trgts.h"
#include "../../DNAQC_trgt_gen/mk_DNA_qc_testdata.h"

#include "../MLT/prcs_MLT_cotrans.h"
#include "../MUX/get_brcd_str.h"

#include "../../DNAQC_trgt_gen/mk_DNA_qc_trgts.h"

#include "./print_DNA_qc_output.h"

//TODO: move to testdata_analysis file once generated
typedef struct testdata_DNA_qc_vars { //structure containing variables for testdata metrics
    int run;                       //flag indicating whether to perform testdata analysis
} testdata_DNA_qc_vars;


/* prcs_rds_DNA_qc: manages processing of TECprobe-MUX data */
int prcs_rds_DNA_qc(TPROBE_names * nm, FILE * fp_MUXtrgs, int trgt_ftype, fastp_params fastp_prms, testdata_DNA_qc_vars * testdata_DNA_qc, int run_bypass_fastp);

/* assess_brcd_concordance: determine whether barcodes present in a sequencing read are concordant */
void assess_brcd_concordance(FILE **ifp, compact_h_node **htbl_MUX, TPROBE_names * nm, compact_target * ctrg, int brcd_cnt, int ctrg_cnt, int * pairing, mapping_metrics * met, int mode);

/* assess_testdata_mapping: assess whether testdata mapping matches expectations */
void assess_testdata_mapping(compact_target * ctrg, int brcd_cnt); 
#endif /* prcs_rds_DNA_qc_h */

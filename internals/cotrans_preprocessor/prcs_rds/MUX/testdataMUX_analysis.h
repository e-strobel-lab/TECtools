//
//  testdataMUX_analysis.h
//  
//
//  Created by Eric Strobel on 7/10/25.
//

#ifndef testdataMUX_analysis_h
#define testdataMUX_analysis_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/gen_utils.h"
#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../MUX_trgt_gen/mk_MUX_trgts.h"

#include "../../../variant_maker/make_barcodes.h"

#define MAX_UINT64_LEN 20 //maximum number of characters in an unsigned 64-bit integer

typedef struct testdata_MUX_vars { //structure containing variables for testdata metrics
    int run;                       //flag indicating whether to perform testdata analysis
    int bid_match;
} testdata_MUX_vars;

/* get_testdata_barcode_id: get barcode id from test data read line 1 */
uint64_t get_testdata_barcode_id(char * id);

/* compare_testdata_ barcode id: check whether test data read mapped to the expected barcode */
void compare_testdata_barcode_id(compact_target * ctrg, char * rd_id, char * sq);

/* print_MUX_testdata_analysis: report outcome of test data read mapping */
void print_MUX_testdata_analysis(metrics * met, compact_target * ctrg);

#endif /* testdataMUX_analysis_h */

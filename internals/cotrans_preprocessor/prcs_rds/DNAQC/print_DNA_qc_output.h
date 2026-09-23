//
//  print_DNA_qc_output.h
//  
//
//  Created by Eric Strobel on 9/16/26.
//

#ifndef print_DNA_qc_output_h
#define print_DNA_qc_output_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/seq2bin_long.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "../../../seq_utils/mapping_metrics.h"

#include "../../../variant_maker/bc_ind.h"

/* print_DNA_qc_output: print DNA QC output files */
void print_DNA_qc_output(char * out_prfx, compact_target * ctrg, int brcd_cnt, int * pairing);

/* merge_bc_id: merge the ids of two barcode targets */
void merge_bc_id(char ** mrg, char * str1, char * str2);

#endif /* print_DNA_qc_output_h */

//
//  mk_BC_TDSPLY_test_data.h
//  
//
//  Created by Eric Strobel on 10/16/25.
//

#ifndef mk_BC_TDSPLY_test_data_h
#define mk_BC_TDSPLY_test_data_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"
#include "../../seq_utils/revcomp.h"

#include "../../variant_maker/variant_maker_defs.h"

#include "../TECdisplay_mapper_structs.h"
#include "../TECdisplay_mapper_defs.h"

#include "../../utils/debug.h"

#include "./mk_TDSPLY_test_data.h"

/* mk_BC_TDSPLY_test_data: coordinates test data generation */
int mk_BC_TDSPLY_test_data(TDSPLY_names * nm, compact_target *ctrg, target_params *trg_prms, int ctrg_cnt);

/* print_BC_TDSPLY_fq: construct read sequences and print to fastq file */
void print_BC_TDSPLY_fq(FILE * out_rd1, FILE * out_rd2, compact_target * ctrg, char * chnl_bc, int end_rnd_typ);

#endif /* mk_BC_TDSPLY_test_data_h */

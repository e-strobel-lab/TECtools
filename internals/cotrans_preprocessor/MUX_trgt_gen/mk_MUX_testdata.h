//
//  mk_MUX_testdata.h
//  
//
//  Created by Eric Strobel on 6/27/25.
//

#ifndef mk_MUX_testdata_h
#define mk_MUX_testdata_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../variant_maker/constant_seqs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../MLT_trgt_gen/mk_3pEnd_testdata.h"

#include "./mk_MUX_trgts.h"

/* mk_MUX_testdata: generate TECprobe-MUX test data */
void mk_MUX_testdata(TPROBE_names * nm, compact_target * ctrg, int ctrg_cnt);

/* print_MUX_fq: print TECprobe-MUX test data read to output fastq file */
int print_MUX_fq(FILE * out_rd1, FILE * out_rd2, char * chnl_bc, compact_target * ctrg);

#endif /* mk_MUX_testdata_h */

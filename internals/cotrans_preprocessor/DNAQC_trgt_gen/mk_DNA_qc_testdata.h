//
//  mk_DNA_qc_testdata.h
//  
//
//  Created by Eric Strobel on 9/10/26.
//

#ifndef mk_DNA_qc_testdata_h
#define mk_DNA_qc_testdata_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../seq_utils/revcomp.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"


/* mk_DNA_qc_testdata: generate TECprobe-MUX test data */
void mk_DNA_qc_testdata(TPROBE_names * nm, compact_target * ctrg, int ctrg_cnt);

/* print_DNA_qc_fq: print TECprobe-MUX test data read to output fastq file */
int print_DNA_qc_fq(FILE * out_rd1, FILE * out_rd2, compact_target * ctrg, compact_target * root_ctrg, int ctrg_cnt, int concordant);

/* find_discordant_ctrg: find barcode ctrg that is discordant relative to the input ctrg*/
compact_target * find_discordant_ctrg(compact_target * root_ctrg, compact_target * ipt_ctrg, int ctrg_cnt);

#endif /* mk_DNA_qc_testdata_h */

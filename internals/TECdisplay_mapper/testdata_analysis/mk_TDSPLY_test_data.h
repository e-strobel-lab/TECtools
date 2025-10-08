//
//  mk_TDSPLY_test_data.h
//  
//
//  Created by Eric Strobel on 4/26/23.
//

#ifndef mk_TDSPLY_test_data_h
#define mk_TDSPLY_test_data_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../variant_maker/variant_maker_defs.h"

#include "../TECdisplay_mapper_structs.h"
#include "../TECdisplay_mapper_defs.h"

#include "../map_reads/map_expected/parse_vmt_trgts.h"

#include "../../utils/debug.h"

/* mk_TDSPLY_test_data: coordinates test data generation */
int mk_TDSPLY_test_data(TDSPLY_names * nm, target *refs, target *trgts, target_params *trg_prms);

/* mk_rndmzd_bc: generate a randomized channel barcode with variable channel and match settings. */
void mk_rndmzd_bc(char * bc, int chnl, int mtch);

/* print_TDSPLY_fq: construct read sequences and print to fastq file */
void print_TDSPLY_fq(FILE * out_rd1, FILE * out_rd2, char * var_id, char * insrt2use, char * chnl_bc, int end_rnd_typ);

/* mutate_insrt: randomly mutate test data sequencing read insert */
int mutate_insrt(char * vb_map, char * mut);

#endif /* mk_TDSPLY_test_data_h */

//
//  cmpr_smlr_seqs.h
//  
//
//  Created by Eric Strobel on 6/29/26.
//

#ifndef cmpr_smlr_seqs_h
#define cmpr_smlr_seqs_h

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"
#include "../../seq_utils/basemap.h"

#include "../../TECdisplay_mapper/TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper/TECdisplay_mapper_structs.h"

#include "../../TECdisplay_mapper/map_reads/map_standard_TDSPLY_reads.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/get_key.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/set_trgt.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"
#include "../../TECdisplay_mapper/map_reads/map_expected/print_navigator_template.h"
#include "../../TECdisplay_navigator/parse_reference.h"
#include "../../cotrans_preprocessor/MUX_trgt_gen/set_barcoded_compact_target.h"
#include "../../variant_maker/constant_seqs.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../set_seq_attributes/parse_descriptor_file.h"

/* opt_sq_att: optional values struct for associating sequence attributes and basemap with a compact_target structure */
typedef struct opt_sq_att {
    sequence_attributes * sq_att;
    basemap * bmap;
} opt_sq_att;


/* cmpr_smlr_seqs: perform a systematic comparison of closely-related sequences to identify sequence variants that exhibit large effects on function */
void cmpr_smlr_seqs(char * msa_path, int msa_typ, sequence_attributes * sq_att, int seq_cnt);

/*set_bmap_from_prsd_vmt: use information from parsed vmt file to generate basemap */
//TODO: need to add ability to include user-specified pairing constraints
void set_bmap_from_prsd_vmt(basemap * bmap, wt_source * wt_src, char ** cnstnt_indels, target * refs, TDSPLY_fasta * wt, target_params * trg_prms);

/* set_sq_att_utl: associate opt_sq_att struct with compact_target struct via the utility member and set values within opt_sq_att struct */
void set_sq_att_utl(compact_target * ctrg, opt_sq_att * utl, sequence_attributes * sq_att, basemap * bmap,  target_params * trg_prms);
#endif /* cmpr_smlr_seqs_h */

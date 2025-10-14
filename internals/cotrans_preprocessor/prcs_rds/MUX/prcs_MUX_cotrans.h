//
//  prcs_MUX_cotrans.h
//  
//
//  Created by Eric Strobel on 6/25/25.
//

#ifndef prcs_MUX_cotrans_h
#define prcs_MUX_cotrans_h

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

#include "../UNV/call_fastp_TPROBE.h"
#include "../UNV/bypass_fastp.h"
#include "../UNV/prcs_chnl_TPROBE.h"
#include "../UNV/print_splitting_metrics.h"
#include "../../../seq_utils/mapping_metrics.h"

#include "../../MUX_trgt_gen/mk_MUX_trgts.h"
#include "../../MUX_trgt_gen/mk_MUX_testdata.h"

#include "../MLT/prcs_MLT_cotrans.h"

#include "./testdataMUX_analysis.h"

#include "prcs_MUX_cotrans.h"

/* prcs_MUX_cotrans: manages processing of TECprobe-MUX data */
int prcs_MUX_cotrans(TPROBE_names * nm, FILE * fp_MUXtrgs, int trgt_ftype, fastp_params fastp_prms, testdata_MUX_vars * testdata_MUX, int run_bypass_fastp);

/* mk_htbl_MUX: makes compact target hash table */
/* hash table has linked list buckets for possible collisions */
void mk_htbl_MUX(compact_h_node ** htbl_MUX, compact_h_node_bank * bank, compact_target * ctrg, int count, mapping_metrics * met);

/* hash_brcd_trgt: generates hash key for binary encoded sequence (up to 32 nt/64 bits */
uint64_t hash_brcd_trgt(binary_seq * bsq);

/* check_brcd_diff: checks whether redundant barcode targets are linked to different barcode variants*/
int check_brcd_diff(compact_target * old, compact_target * new);

/* map_brcd: map barcode to target using hash table */
compact_target * map_brcd(char * brcd_str, char * rc_brcd_str, compact_h_node **htbl_MUX, compact_target ** mpd_trg, mapping_metrics * met);

/* split_MUX_reads: demultiplex TECprobe-MUX reads into separate fastq file */
void split_MUX_reads(FILE **ifp, compact_h_node **htbl_MUX, TPROBE_names * nm, compact_target * ctrg, int brcd_cnt, int ctrg_cnt, mapping_metrics * met, int mode);

/* get_brcd_str: get barcode string from UMI in read ID */
void get_brcd_str(char * brcd_str, char * read1_ID);

#endif /* prcs_MUX_cotrans_h */

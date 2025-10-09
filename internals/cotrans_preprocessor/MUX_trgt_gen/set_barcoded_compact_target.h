//
//  set_barcoded_compact_target.h
//  
//
//  Created by Eric Strobel on 10/1/25.
//

#ifndef set_barcoded_compact_target_h
#define set_barcoded_compact_target_h

#include <stdio.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"

#include "../../TECdisplay_mapper/map_reads/map_expected/parse_vmt_trgts.h"

#define NUMERICAL_ID 0    //indicates that barcode ID is strictly numerical
#define COMPLEX_ID 1      //indicates that barcode ID is not strictly numerical

#define NATIVE_BRCD 0     //barcode contains no mutations
#define MUTANT_BRCD 1     //barcode contains a mutation

#define MUTCODE_BITS 16   //number of bits allocated for barcode mutation codes

/* set_barcoded_compact_target: set target values in target struct */
void set_barcoded_compact_target(compact_target * ctrg, opt_BC * BC_val, target * crnt_ref, char * trgt_id, char * trgt_sq, target_params * trg_prms, int trgt_ftype);

/* parse_barcode_id: parse barcode id to set full id and numerical id pointers */
int parse_barcode_id(char ** p_nid, char ** p_fid, char * crnt_bcid, int trgt_ftype);

/* parse_target_seq: parse oligonucleotide sequence to identify barcode and target RNA sequences */
void parse_target_seq(char ** p_trgt, char ** p_brcd, int * brcd_len, char * crnt_seq, int trgt_ftype);

/* set_BC_val: set optional values for barcode target */
void set_BC_val(compact_target * ctrg, opt_BC * BC_val, char * tsq, compact_target * ntv, int mode);

#endif /* set_barcoded_compact_target_h */

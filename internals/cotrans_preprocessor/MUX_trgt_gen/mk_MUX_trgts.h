//
//  mk_MUX_trgts.h
//  
//
//  Created by Eric Strobel on 6/25/25.
//

#ifndef mk_MUX_trgts_h
#define mk_MUX_trgts_h

#include <stdio.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/constant_seqs.h"
#include "../../variant_maker/read_bcFile.h"

#define MAX_UINT64_LEN 20 //maximum number of characters in an unsigned 64-bit integer

#define MUTCODE_BITS 16   //number of bits allocated for barcode mutation codes
#define NUMERICAL_ID 0    //indicates that barcode ID is strictly numerical
#define COMPLEX_ID 1      //indicates that barcode ID is not strictly numerical

#define NATIVE_BRCD 0     //barcode contains no mutations
#define MUTANT_BRCD 1     //barcode contains a mutation

#define MIN_SUB_CODE 1    //minimum mutcode value for substitutions
#define MAX_SUB_CODE 48   //maximum mutcode value for substitutions
#define MIN_INS_CODE 49   //minimum mutcode value for insertions
#define MAX_INS_CODE 112  //maximum mutcode value for insertions
#define MIN_DEL_CODE 113  //minimum mutcode value for deletions
#define MAX_DEL_CODE 128  //maximum mutcode value for deletions

/* mk_MUX_trgts: manages TECprobe-MUX target generation */
int mk_MUX_trgts(compact_target * ctrg, opt_BC * BC_val, FILE * fp_MUXtrgs, int brcd_cnt, int clcd_ctrg_cnt);

/* store_barcode_targets: stores input barcodes as compact targets */
int store_barcode_targets(compact_target * ctrg, opt_BC * BC_val, int * brcd_len, FILE * fp_MUXtrgs, int brcd_cnt);

/* set_BC_val: set optional values for barcode target */
void set_BC_val(compact_target * ctrg, opt_BC * BC_val, char * tsq, compact_target * ntv, int mode);

/* parsce_barcode_id: parse barcode id to set full id and numerical id pointers */
int parse_barcode_id(char ** p_nid, char ** p_fid, char * crnt_bcid);

/* parse_oligo_seq: parse oligonucleotide sequence to identify barcode and target RNA sequences */
void parse_oligo_seq(char ** p_trgt, char ** p_brcd, int * brcd_len, char * crnt_oligo);

/* mk_SUB_trgts: generates single substitution targets for input barcode */
int mk_SUB_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode);

/* mk_INS_trgts: generate single insertion targets for input barcode */
int mk_INS_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode);

/* mk_DEL_trgts: generate single deletion targets for input barcode */
int mk_DEL_trgts(compact_target * ctrg, opt_BC * BC_val, int * ctrg_cnt, compact_target * src_ctrg, int brcd_len, uint64_t * mutCode, char upstrm_nt);

/* get_target_type: determine target type using mutcode */
char * get_target_type(uint64_t mutcode, int * type_val);

#endif /* mk_MUX_trgts_h */

//
//  parse_DNA_qc_trgts.h
//  
//
//  Created by Eric Strobel on 9/8/26.
//

#ifndef parse_DNA_qc_trgts_h
#define parse_DNA_qc_trgts_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../../variant_maker/make_barcodes.h"
#include "../../variant_maker/constant_seqs.h"
#include "../../variant_maker/bc_ind.h"

#include "../MUX_trgt_gen/set_barcoded_compact_target.h"

#include "./link_trgts.h"

/* parse_DNA_qc_trgts: parse DNA qc targets fasta file */
void parse_DNA_qc_trgts(FILE * ifp, int trgt_ftype, compact_target * ctrg, opt_BC * trg_val, association * assoc, target_params * trg_prms);

/* parse_DNA_qc_trgt_id: parse target id of DNA qc target */
void parse_DNA_qc_trgt_id(char ** var_id, char ** bc_id1, char ** bc_id2, char * trgt_id);

/* parse_DNA_qc_barcodes: parse barcode sequences of DNA qc target */
int parse_DNA_qc_barcodes(char ** bc1, char ** bc2, char * prsd_bc_sq);

/* parse_DNA_qc_insert: parse insert of DNA qc target */
void parse_DNA_qc_insert(char **insert, char * prsd_insert_sq, int brcd_len);

/* set_DNA_qc_compact_target: set compact_target values for DNA qc target */
void set_DNA_qc_compact_target(compact_target * ctrg, opt_BC * BC_val, char * insert, char * var_id, char * bc_id, char * bc_sq, int bc_num, target_params * trg_prms);

#endif /* parse_DNA_qc_trgts_h */

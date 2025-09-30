//
//  parse_mx_trgts.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef parse_mx_trgts_h
#define parse_mx_trgts_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/revcomp.h"
#include "../../../seq_utils/isDNAbase.h"
#include "../../../seq_utils/isIUPACbase.h"
#include "../../../seq_utils/seq2bin_hash.h"
#include "../../../seq_utils/basemap.h"

#include "../map_reads.h"

#define TDSPLY_TRGS 0
#define TPROBE_TRGS 1

/* parse_header_lines: parse targets file header lines for expected variant count
 and wild type sequence information */
void parse_header_lines(FILE * ifp, target_params *trg_prms, fasta * wt);

/* parse_mx_trgts: parse targets file to obtain target ids, sequences, attributes,
  and min/max transcript lengths */
void parse_mx_trgts(FILE * ifp, target * refs, opt_ref * ref_val, void * trgts, void * optmx, target_params * trg_prms, fasta * wt, int mode);

/* check_tpr_match: check that expected tpr matches actual tpr */
void check_tpr_match(int cnt, int actual, int xpctd);

/* validate_ref_seq: check that reference comprises valid characters and is not too long */
void validate_ref_seq(char * ref_nm, char * ref_sq, fasta * wt);

/* validate_trgt_seq: check that target comprises valid characters and
 matches its associated reference target at all non-variable positions */
int validate_trgt_seq(char * trgt_nm, char * trgt_sq, char * ref_sq);

/* process_trgt_seq: convert target seq to upper case and remove non-base characters. */
int process_trgt_seq(char * ipt, char * processed_seq);

/* set_ref: set reference target values in targets structure */
void set_ref(target * refs,  opt_ref * ref_val, char * ref_id, char * ref_sq, int tpr, char * vbases, char * cnstnts);

/* set_trgt: set target values in target struct */
void set_trgt(target * trgts, opt_mx_trg * trg_val, target * crnt_ref, char * trgt_id, char * trgt_sq);

/* print_reference_debug: print debug messages for reference processing outcome */
void print_reference_debug(target * trgt, target_params * trg_prms);

/* print_target_debug: print debug messages for target processing outcome */
void print_target_debug(target * trgt, target_params * trg_prms);

#endif /* parse_mx_trgts_h */

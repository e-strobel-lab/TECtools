//
//  prcs_MLT_cotrans.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef prcs_MLT_cotrans_h
#define prcs_MLT_cotrans_h

#include <stdio.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/seq2bin_hash.h"
#include "parse_3pEnd_trgts.h"
#include "../UNV/call_fastp.h"
#include "../UNV/prcs_chnl.h"
#include "../../../seq_utils/appnd_att.h"
#include "../../../seq_utils/isDNAbase.h"
#include "printQC_prcsMLT.h"
#include "testdata3pEnd_analysis.h"
#include "mk_smooth_script.h"
#include "mk_MLT_config.h"
#include "../../../utils/debug.h"


/* prcs_MLT_cotrans: prcs_MLT_cotrans manage the functions that are required to split fastq files from
 multi-length cotranscriptional RNA structure probing experiments into separate fastq files based on the
 channel barcode and the RNA 3' end.
 
 ***arguments***
 names * nm: contains input R1 and R2 fastq filenames and 3' end targets filename
 FILE * fp_3pEND: pointer the 3' end targets file
 fastp_params fastp_prms: parameters for fastp processing
 
 */
int prcs_MLT_cotrans(names * nm, FILE * fp_3pEnd, fastp_params fastp_prms);



/* mk_htbl_3pEnd: generate a hash table for 3' end targets using seq2bin_hash functions
 
 ***arguments***
 h_node **htbl: pointer to hash table node pointer array
 h_node_bank *bank: pointer to bank of hash table nodes
 target *trgts: pointer to array of hash table target structures
 int count: the number of 3' end targets
 */
void mk_htbl_3pEnd(h_node **htbl, h_node_bank *bank, target *trgts, int count); //construct hash table



/* check_diff: assesses hash table collisions by checking the difference in transcript length
 between the previous hash table target and the new target is less than TRG_DIFF_THRESHOLD
 (defined in cotrans_preprocessor_structs.h)
 
 ***arguments***
 target prev: target that was already entered into the hash table
 target new: target with a sequence that matches that of prev
 */

void check_diff(target prev, target new);


/* map_3pEnd: use the 3' end targets hash table to determine the transcript length of a read
 
 ***arguments***
 char * read: full sequence of the read that will be assessed (always read 1 of a pair)
 h_node **htbl: pointer to hash table node pointer array
 char * end_str: array for storing transcript length value as chars
 metrics * met: pointer metrics structure for read processing QC
 int trg_len: length of 3' end target sequences
 */
int map_3pEnd(char * read, h_node **htbl, char * end_str, metrics * met, int trg_len);


/* verify_read: perform checks to assess the integrity of a read */
int verify_read(char (*rd)[MAX_LINE]);

/* split_reads_3pEnd: split input fastq files by channel barcode and 3' end.
 
 ***arguments***
 FILE **ifp: pointer to read 1 and read 2 fastq file pointers
 h_node **htbl:  pointer to hash table node pointer array
 names * nm: contains sample names to use as prefix when generating output fastq files
 metrics  * met: pointer metrics structure for read processing QC
 target3p_params trg_prms: 3' end targets parameters
 int mode: MULTI or SINGLE mode specification
 */

int split_reads_3pEnd(FILE **ifp, h_node **htbl, names * nm, metrics  * met, target3p_params trg_prms, int mode);


#endif /* prcs_MLT_cotrans_h */

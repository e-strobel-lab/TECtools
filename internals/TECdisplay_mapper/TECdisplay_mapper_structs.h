//
//  TECDisplay_mapper_structs.h
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#ifndef TECdisplay_mapper_structs_h
#define TECdisplay_mapper_structs_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/fastp_params.h"

#include "./TECdisplay_mapper_defs.h"
#include "../seq_utils/seq2bin_hash.h"

#include "../variant_maker/variant_maker_defs.h"

/* NOTE: the structures used for the hash table are defined in seq2bin_hash.h */

/* opt_ref: structure for storing reference target values.
 this structure is used as the 'opt' member of the target structure defined
 in seq2bin_hash.h for reference targets */
typedef struct opt_ref {
    char * ipt_sq;                 //input sequence (contains spacers) for printing to navigator template
    char * vbases;                 //variable base string for constructing navigator template
    char * cnstnts;                //constant insertions and deletions string
    int vb_pos[SEQ2BIN_MAX_KEY+1]; //position (indices) of variable bases in the reference sequence
    int vb_cnt;                    //number of variable bases in the reference sequence
    int tpr;                       //number of targets encoded in ref
} opt_ref;


/* opt_mx_trg: structure for storing TECdisplay target values.
 this structure is used as the 'opt' member of the target structure defined
 in seq2bin_hash.h for target sequences */
typedef struct opt_mx_trg {
    struct target * ref;            //pointer to reference sequence for the target
    char vb_key[SEQ2BIN_MAX_KEY+1]; //hash key string composed of variable bases within the target
    int bnd;                        //number of reads that mapped to the target with a BOUND barcode
    int unb;                        //nubmer of reads that mapped to the target with an UNBOUND barcode
} opt_mx_trg;


/* target_params: structure for storing parameters of TECdisplay targets.
 values in the structure are populated by the parse_mx_trgts function */
typedef struct target_params {
    int xpctd;                      //expected number of targets based on targets file first line
    int r_cnt;                      //number of reference targets
    int t_cnt;                      //number of targets
    int t_per_r[MAXREF];            //number of targets for each reference target
    int nr_cnt;                     //number of non-redundant targets
    int rdndnt;                     //number of redundant targets
    int mapped2;                    //number of targets with at least 1 mapped read
    int BClen;                      //barcode length (used in some cases)
} target_params;


/* names: structure containing key names used during seq read processing */
typedef struct TDSPLY_names {
    char file[READ_MAX][MAX_LINE];  //array to store read 1 and 2 filenames
    char smpl[READ_MAX][MAX_LINE];  //array to store read 1 and 2 sample names
    char mrg[MAX_LINE];             //array to store merged read sample name
    char trgs[MAX_LINE];            //array to store targets filename
    char trgs_prefix[MAX_LINE];     //array to store targets prefix
    char out_nm[MAX_LINE];          //array to store output file name
} TDSPLY_names;

/* fasta: structure to store fasta entries */
typedef struct TDSPLY_fasta {
    char * nm;           //sequence name
    char * sq;           //sequence with multiple sequence alignment characters preserved
    int * pos;           //position number of base
    int ex;              //flag to exclude variant
    int start_pos;       //position of first nucleotide in sequence
} TDSPLY_fasta;

/* metrics: structure containing seq read processing metrics */
typedef struct TDSPLY_metrics {
    int reads_processed;            //number of processed reads
    int chan_count[CHANNEL_MAX];    //tracks number of reads from each channel
    int full_match[CHANNEL_MAX];    //tracks number of full channel matches
    int part_match[CHANNEL_MAX];    //tracks number of partial channel matches
    int hits;                       //tracks number of reads wit ha key that maps to the hash table
    int matches;                    //tracks number of reads with an exact target match
} TDSPLY_metrics;


/* testdata_vars: structure containing variable for tracking test data processing metrics */
typedef struct testdata_vars {      //structure containing variables for testdata metrics
    int run;                        //flag indicating whether to perform testdata analysis
    int nat_rds;                    //number of native reads in test data input
    int nat_mpd;                    //number of native reads that mapped to a target
    int nat_exp;                    //number of native reads that mapped to the expected target
    int nat_err;                    //number of native reads that mapped to an unexpeced target
    int mut_rds;                    //number of mutant reads in test data input
    int mut_mpd;                    //number of mutant reads that mapped to a target
} testdata_vars;

#endif /* TECdisplay_mapper_structs_h */

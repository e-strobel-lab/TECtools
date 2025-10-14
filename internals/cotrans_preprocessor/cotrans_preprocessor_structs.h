//
//  cotrans_preprocessor_structs.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef cotrans_preprocessor_structs_h
#define cotrans_preprocessor_structs_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/fastp_params.h"
#include "../seq_utils/seq2bin_hash.h"
#include "../seq_utils/seq2bin_long.h"

#include "../TECdisplay_mapper/TECdisplay_mapper_structs.h"

#include "cotrans_preprocessor_defs.h"

/* NOTE: the structures used for the hash table are defined in seq2bin_hash.h */

/* target3p_genVals: structure containing variables for end target generation */
typedef struct target3p_genVals {
    char sq[512][MAX_END_LEN+2];	//array for native end target sequences
    int len;						//length of end target sequences
    int nat;						//number of native end variants
    int dstnc[MAX_END_LEN+1];		//array to track distribution of end target ambiguity
    int exp[4];						//expected number of NAT, SUB, INS, and DEL end targets
    int act[4];						//actual number of NAT, SUB, INS, and DEL end targets
} target3p_genVals;


/* opt_3pEnd: structure for storing the length and type of a 3' end target.
 this structure is used as the 'opt' member of the target structure defined
 in seq2bin_hash.h when processing multi-length cotrans data */
typedef struct opt_3pEnd {
    int len;						//transcript length that the 3' end target corresponds to
    int typ;						//code for type of target {NAT, SUB, INS, or DEL}
} opt_3pEnd;


/* target3p_params: structure for storing parameters of 3' end targets.
 values in the structure are populated by the parse_3pEnd_trgts function */
typedef struct target3p_params {
    int min;						//minimum transcript length
    int max;						//maximum transcript length
    int sdLen;						//length of seed sequence for 3' end mapping
    int cnt;						//number of targets
} target3p_params;


/* opt_BC: structure for storing barcode target-specific variables */
typedef struct opt_BC {
    compact_target * ref;                     //pointer to reference target, for native targets this points to self
    FILE * ofp[TPROBE_CHANNEL_MAX][READ_MAX]; //output file pointer (only used for native targets)
    uint64_t mpd;                      //number of reads that mapped to native+mutant targets (used for native targets)
    uint64_t cnt;                      //number of reads that mapped to the current target
    char * tsq;                        //template sequence
    int typ;                           //target type (NAT, SUB, INS, DEL)
} opt_BC;


/* TPROBE_names: structure containing key names during seq read processing */
typedef struct TPROBE_names {
    char file[READ_MAX][MAX_LINE+1]; //array to store read 1 and 2 filenames
    char smpl[READ_MAX][MAX_LINE+1]; //array to store read 1 and 2 sample names
    char trgts[MAX_LINE+1];          //array to store 3' end targets filename
    char trgts_prfx[MAX_LINE+1];     //barcode filename prefix
    int len;                         //length of single length target
} TPROBE_names;

#endif /* cotrans_preprocessor_structs_h */

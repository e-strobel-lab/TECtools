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


/* names: structure containing key names during seq read processing */
typedef struct names {
    char file[READ_MAX][MAX_LINE];	//array to store read 1 and 2 filenames
    char smpl[READ_MAX][MAX_LINE];	//array to store read 1 and 2 sample names
    char ends[MAX_LINE];			//array to store 3' end targets filename
    int len;                        //length of single length target
} names;


/* metrics: structure containing seq read processing metrics */
typedef struct metrics {
    int reads_processed;			//number of processed reads
    int chan_count[CHANNEL_MAX];	//count of reads from each channel
    int full_match[CHANNEL_MAX];	//counts number of full channel matches
    int part_match[CHANNEL_MAX];	//counts number of partial channel matches
    int end_hits;					//tracks number of ends that were mapped by hash table
    int end_matches;				//tracks number of end hits with expected sequence (should be 100%)
    int native_cnt;					//tracks number of observed 3' ends that match native targets
    int sub_cnt;					//tracks number of observed 3' ends that match substitution targets
    int mapped_ends;				//mapped ends count
    int unmapped_ends;				//count of ends that failed to map
    int read_matches[READ_MAX];		//count of post-processing verified reads
    int len_dist[MAX_LINE];			//array for tracking the length distribution of UNT, MOD, and ERR reads
} metrics;


/* fastp_params: structure containing parameters for fastp processing */
typedef struct fastp_params {
    char path[MAX_LINE];            //path to fastp executable
    int mode;						//processing mode (multi-length or single-length)
    int limit;						//limit on number of reads to process
} fastp_params;						//NOTE: processing mode is also used by functions unrelated to fastp



#endif /* cotrans_preprocessor_structs_h */

//
//  cotrans_preprocessor_defs.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef cotrans_preprocessor_defs_h
#define cotrans_preprocessor_defs_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../TECdisplay_mapper/TECdisplay_mapper_defs.h"

/* run_modes */
#define MK_3pEND_TARGETS 0	//multi-length cotrans target generation mode
#define PRCS_READS 1		//seq read pre-processing mode
#define MK_RUN_SCRIPT 2		//shapemapper2 run script generation mode
#define MK_FASTA 3          //fasta generation mode
#define MK_SINGLE_TARGET 4  //single-length cotrans target generation mode

/* 3pEnd targets definitions */
#define DFLT_MIN 19			//min length is 19 to include seq information for length 20 deletions
#define END_LEN_INIT 14 	//default end target length
#define MAX_END_LEN 20		//maximum end target length, must be <MAX_TRGT_SQ
#define MATCHDIF 3			//threshold for assessing end target ambiguity
#define HEADER 0			//code for printing header with print_end_target_table function
#define LINE 1				//code for printing data line from print_end_target_table function todo rename
#define FULL 0				//full match code
#define PART 1				//partial match code
#define CHNL_CLASSES 5		//number of channel classes to generate when making test data

/* read processing modes */
#define MULTI 0				//multi-length cotranscriptional mode
#define SINGLE 1			//single length mode
#define MULTIPLEX 2         //multiplex mode

/* target parsing definitions */
#define CNSTNT_STRT_NAT 7   //offset from last target id string to reach start of standardized target id (NAT)
#define CNSTNT_STRT_MUT 11  //offset from last target id string to reach start of standardized target id (MUT)

/* end mapping definitions */
#define END_MAX 999			//maximum allowable 3' end length
#define MAX_TRGT_ID 32		//maximum target id length
#define MAX_TRGT_SQ 32		//maximum target sequence length, must be >MAX_END_LEN
#define LEN_CODE 4	//array size for storing length info when parsing 3' end files
#define MUT_CODE 4	//array size for storing match type information when parsing 3' end files
#define SUB 0				//code for match to 3' end sequence with 1 substitution
#define INS 1				//code for match to 3' end sequence with 1 insertion
#define DEL 2				//code for match to 3' end sequence with 1 deletion
#define NAT 3				//code for match to native 3' end sequence
#define TRG_DIF_THRESHOLD 2	//max distance for 3' value of targets with the same sequence

/* channel ID definitions */
#define VL_CHNL_BC_LEN 5	//length of the channel barcode at the head of read 1
#define VL_MIN_MATCH 4	    //minimum channel barcode match (4/5 nucleotides)
#define VL_MAX_MATCH 5	    //maximum possible channel barcode match (5/5 nucleotides)
#define VL_UMI_LENGTH 9		//length of the unique molecular identifier
#define UNT 0				//specifies untreated sample
#define MOD 1				//specifies modified sample
#define ERR 2				//specifies undetermined channel barcode
#define CHANNEL_MAX 3		//number of channel codes (UNT, MOD, ERR)

/* barcode mapping definitions */
#define MAX_BRCD_CNT 20000  //maximum number of TECprobe-MUX barcodes

/* split output definitions */
#define FILES_PER_LEN 4		//number of files generated per transcript length

#endif /* cotrans_preprocessor_defs_h */

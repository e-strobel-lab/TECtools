//
//  TECDisplay_mapper_main.h
//  
//
//  Created by Eric Strobel on 6/21/22.
//

#ifndef TECdisplay_mapper_main_h
#define TECdisplay_mapper_main_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../variant_maker/variant_maker_defs.h"
#include "../seq_utils/seq2bin_hash.h"

/* run_modes */
#define MAP_TEST_DATA 0      //map test data mode
#define MAP_SEQ_READS 1      //map sequencing reads mode

/* quality score threshold indices */
#define Q_VARIABLE 0         //index of minimum variable base quality score
#define Q_CONSTANT 1         //index of minimum constant base quality score

/* channel ID definitions */
#define CHANNEL_BC_LENGTH 4  //length of the channel barcode at the head of read 2
#define MIN_MATCH 3          //minimum channel barcode match (3/4 nucleotides)
#define MAX_MATCH 4          //maximum possible channel barcode match (4/4 nucleotides)
#define UMI_LENGTH 16        //length of the unique molecular identifier
#define BND 0                //specifies bound sample
#define UNB 1                //specifies unbound sample
#define ERR 2                //specifies undetermined channel barcode
#define CHANNEL_MAX 3        //number of channel codes (BND, UNB, ERR)

/* test data generation */
#define NO_CHNL_MATCH -1     //no match code
#define FULL 0               //full match code
#define PART 1               //partial match code
#define CHNL_CLASSES 5       //number of channel classes to generate when making test data

#define SUB 0                //code for sequence with 1 substitution
#define INS 1                //code for sequence with 1 insertion
#define DEL 2                //code for sequence with 1 deletion
#define NAT 3                //code for native 3' sequence

#define READ_KEY 0           //flag that input sequence for get_key comes from a read
#define TARGET_KEY 1         //flag that input sequence for get_key comes from a target

#define REFERENCE 0          //flag to print target processing debug for reference target
#define TARGET 1             //flag to print target processing debug for target

/* navigator template generation */
#define BASE_FIELD 8         //width of base constraint field in navigator template file


#endif /* TECdisplay_mapper_main_h */

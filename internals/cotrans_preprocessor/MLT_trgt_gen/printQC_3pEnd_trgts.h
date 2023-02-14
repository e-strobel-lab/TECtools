//
//  printQC_3pEnd_trgts.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef printQC_3pEnd_trgts_h
#define printQC_3pEnd_trgts_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/gen_utils.h"

/* print_native_end_QC: prints table that contains QC information for each native end target
 
 ***arguments***
 FILE * metrics_fp: pointer to 3' end metrics output file
 target3p_genVals * ends: end target generation values
 int len: 3' end length of table line
 int fndMtch: number of matches with distance <MATCHDIF found
 int mtchPos: 3' end length of exact match, 99999 if more than one exact match
 int mode: specifies whether to print header or data line
 */
int print_native_end_QC(FILE * metrics_fp, target3p_genVals * ends,
                        int len, int fndMtch, int mtchPos, int mode);



/* print_match_distance_table: prints table that contains summary data
 describing distance of each end target from all other end targets
 
 ***arguments***
 FILE * metrics_fp: pointer to 3' end metrics output file
 target3p_genVals * ends: end target generation values
 int total_tested: total number of native sequences tested
 */
int print_match_distance_table(FILE * metrics_fp, target3p_genVals * ends, int total_tested);



/* print_end_target_table: prints table that contains the number of end targets generated
 ***arguments***
 
 FILE * metrics_fp: pointer to 3' end metrics output file
 target3p_genVals * ends: end target generation values
 */
int print_end_target_table(FILE * metrics_fp, target3p_genVals * ends);

#endif /* printQC_3pEnd_trgts_h */

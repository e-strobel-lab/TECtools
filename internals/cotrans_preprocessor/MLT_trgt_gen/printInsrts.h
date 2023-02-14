//
//  printInsrts.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef printInsrts_h
#define printInsrts_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

/* printInsrts: generate all possible 1 nt insertions for each end target
 
 internal insertions can be any of the four DNA bases because there is not
 a native nucleotide at this position (term: (ipt_len - min_len) * 12 * 4).
 leading and trailing insertions can only be the 3 non-native DNA bases
 because inserting the native base at the leading position will generate a
 sequence that appears the same as the native target and inserting the native
 base at the trailing position will appear the same as the next native target
 (term: (ipt_len - min_len) * 6). the exception to this rule is the trailing
 insertion of the last 3' end target (term: + 1).
 
 calculated value: ((ipt_len - min_len) * 6) + ((ipt_len - min_len) * 12 * 4) + 1
 
 ***arguments***
 target3p_genVals * ends: end target generation values
 char name[MAX_LINE]: sequence namee
 FILE * out_fp: pointer to 3' end targets output file
 int min_len: minimum transcript length-1
 */
int printInsrts(target3p_genVals * ends, char name[MAX_LINE], FILE * out_fp, int min_len);

#endif /* printInsrts_h */

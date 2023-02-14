//
//  printSubs.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef printSubs_h
#define printSubs_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

/* printSubs: generate all possible 1 nt substitutions for each end target
 
 starting at length 30, substitute each nucleotide
 with the other three possible nucleotides
 
 calculated value: (ipt_len - min_len) * len * 3
 
 ***arguments***
 target3p_genVals * ends: end target generation values
 char name[MAX_LINE]: sequence namee
 FILE * out_fp: pointer to 3' end targets output file
 int min_len: minimum transcript length-1
 */
int printSubs(target3p_genVals * ends, char name[MAX_LINE], FILE * out_fp, int min_len);

#endif /* printSubs_h */

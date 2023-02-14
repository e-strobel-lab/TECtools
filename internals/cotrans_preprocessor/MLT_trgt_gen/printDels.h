//
//  printDels.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef printDels_h
#define printDels_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

/* printDels: generate all possible 1 nt deletions for each end target
 
 starting at length 30, delete each nucleotide and
 append the first character of the preceding target
 to the start of the current target
 
 it is not possible to delete the last nucleotide in
 a target because this would appear identical to the
 preceding target
 
 calculated value: (ipt_len - min_len) * (len - 1)
 
 ***arguments***
 target3p_genVals * ends: end target generation values
 char name[MAX_LINE]: sequence namee
 FILE * out_fp: pointer to 3' end targets output file
 int min_len: minimum transcript length-1
 */
int printDels(target3p_genVals * ends, char name[MAX_LINE], FILE * out_fp, int min_len);

#endif /* printDels_h */

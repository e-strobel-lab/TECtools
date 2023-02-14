//
//  mk_3pEnd_trgts.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef mk_3pEnd_trgts_h
#define mk_3pEnd_trgts_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "printSubs.h"
#include "printInsrts.h"
#include "printDels.h"
#include "printQC_3pEnd_trgts.h"

/* mk_3pEnd_trgts: generate a file containing every possible 3' end
 of length len, including native sequences, 1 nt substitutions,
 1 nt insertions, and 1 nt deletions
 */
int mk_3pEnd_trgts(target3p_genVals * ends, char name[MAX_LINE], char seq[MAX_LINE], int min_len);

#endif /* mk_3pEnd_trgts_h */

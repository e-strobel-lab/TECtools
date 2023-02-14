//
//  mk_intermed_trgts.h
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#ifndef mk_intermed_trgts_h
#define mk_intermed_trgts_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

/* mk_intermed_trgts: generate individual fasta files for every intermediate
 length of the input sequence, starting at min_len+1 */
int mk_intermed_trgts(char name[MAX_LINE], char seq[MAX_LINE], int min_len);

#endif /* mk_intermed_trgts_h */

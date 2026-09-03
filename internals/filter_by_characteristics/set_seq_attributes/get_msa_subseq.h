//
//  get_msa_subseq.h
//  
//
//  Created by Eric Strobel on 11/24/25.
//

#ifndef get_msa_subseq_h
#define get_msa_subseq_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../seq_utils/isIUPACbase.h"

#define BOUND2_INDEX  0
#define BOUND2_LENGTH 1

/* get_msa_subseq: get subsequence from multiple sequence alignment line */
//NOTE: input sequence indexing must match the true numbering of the sequence (index 0 should have a placeholder char)
int get_msa_subseq(char ** seq, char * msa, int b1, int b2, int mode);

#endif /* get_msa_subseq_h */

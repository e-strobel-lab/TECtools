//
//  get_key.h
//  
//
//  Created by Eric Strobel on 9/30/25.
//

#ifndef get_key_h
#define get_key_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../utils/debug.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"

/* get_key: generate key string composed the nucleotides at variable
 base positions in the input read sequence */
int get_key(char * key, char * end5p, char * qscore5p, char * minQv, target *refs, int key_type);

#endif /* get_key_h */

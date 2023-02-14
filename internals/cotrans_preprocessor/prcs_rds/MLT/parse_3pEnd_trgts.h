//
//  parse_3pEnd_trgts.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef parse_3pEnd_trgts_h
#define parse_3pEnd_trgts_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"
#include "../../../utils/io_management.h"
#include "../../../seq_utils/revcomp.h"
#include "../../../seq_utils/seq2bin_hash.h"

/* parse_3pEnd_trgts: parse targets file to obtain target ids, sequences, attributes,
 and min/max transcript lengths
 
 ***arguments***
 FILE * ifp: pointer to 3' end targets input file
 target * trgts: pointer to array for storing hash table targets
 opt_3pEnd * end3p: pointer to array for storing optional 3' end values;
    becomes linked to the 'opt' member of the corresponding target structure
 target3p_params * trg_prms: 3' end target parameters
 */

void parse_3pEnd_trgts(FILE * ifp, target * trgts, opt_3pEnd * end3p, target3p_params * trg_prms);

#endif /* parse_3pEnd_trgts_h */

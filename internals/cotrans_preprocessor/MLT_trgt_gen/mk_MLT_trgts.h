//
//  mk_MLT_trgts.h
//  
//
//  Created by Eric Strobel on 3/16/22.
//

#ifndef mk_MLT_trgts_h
#define mk_MLT_trgts_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"
#include "../../utils/io_management.h"
#include "../../seq_utils/parse_fasta.h"
#include "mk_intermed_trgts.h"
#include "mk_3pEnd_trgts.h"
#include "mk_3pEnd_testdata.h"

/* mk_MLT_trgts: manage intermediate and 3' end target generration
 ***arguments***
 FILE * fp_fasta: pointer to input fasta file
 char * trgts_nm: pointer to user-supplied target namee
 int usr_min_len: user-supplied minimum target length
 int usr_end_len: user-supplied minimum 3' end target length
 int make_test_data: flag to generate test data
 int mltplir: multiplier value for generating more than the default amount of test data
 int randomize_end: flag to mutate 3' end region when making test data
 */
int mk_MLT_trgts(FILE * fp_fasta, int usr_min_len, int usr_end_len, int make_test_data, int mltplir, int randomize_end);

/* set_end_len: set the length of 3' end target sequences
 int usr_end_len: user-supplied minimum end length
 target3p_genVals * ends: pointer to 3' end generation variables structure
 */
int set_end_len(int usr_end_len, target3p_genVals * ends);

/* set_min_len: the length of the shortest target
 int usr_min_len: user-supplied minimum target length
 int * min_len: pointer minimum target length that will be used for target generation
 */
int set_min_len(int usr_min_len, int * min_len);

#endif /* mk_MLT_trgts_h */

//
//  mk_SGL_target.h
//  
//
//  Created by Eric Strobel on 1/17/23.
//

#ifndef mk_SGL_target_h
#define mk_SGL_target_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "../../utils/io_management.h"
//mk_MLT_targets is only included for the parse_fasta
//function. should move parse_fasta elsewhere
#include "../MLT_trgt_gen/mk_MLT_trgts.h"

/* make fasta file for single length cotrans target */
int mk_SGL_target(FILE * fp_fasta);
 
#endif /* mk_SGL_target_h */

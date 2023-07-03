//
//  read_varFile.h
//  
//
//  Created by Eric Strobel on 8/3/22.
//

#ifndef read_varFile_h
#define read_varFile_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/io_management.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"


/* read_varFile: open file input and assign to file pointer */
int read_varFile(FILE * ifp, wt_source * fa_ipt, basemap * bmap, int mode);

/* validate_seq_line: check that sequence line is valid */
int validate_seq_line(char * id, char * seq, int * id_lkup, int wtsq_len);

/* validate_pair_line: check that pair line is valid */
int validate_pair_line(char * id_p, char * id_s, char * prs, int wtsq_len);

#endif /* read_varFile_h */

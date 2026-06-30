//
//  print_N_compressed_fa.h
//  
//
//  Created by Eric Strobel on 6/19/26.
//

#ifndef print_N_compressed_fa_h
#define print_N_compressed_fa_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

/* print_N_compressed_fa: print fasta file in which variants are compressed by replacing nucleotides that can be any base with 'N' */
void print_N_compressed_fa(names * nm, basemap * bmap, int vTmpCnt, int varCnt, char * out_dir, int append_priming, int append_barcode, char * lnkr, int make_fasta, int lib_type);

#endif /* print_N_compressed_fa_h */

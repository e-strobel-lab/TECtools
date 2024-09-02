//
//  read_bcFile.h
//  
//
//  Created by Eric Strobel on 8/28/24.
//

#ifndef read_bcFile_h
#define read_bcFile_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/io_management.h"
#include "../utils/gen_utils.h"
#include "../seq_utils/isDNAbase.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#define HEADER 0   //parsing barcode file header
#define BC_LINE 1  //parsing barcode line

/* read_bcFile: parse barcode file header and barcode lines */
int read_bcFile(FILE * fp_brcd, int mode, char * str);

#endif /* read_bcFile_h */

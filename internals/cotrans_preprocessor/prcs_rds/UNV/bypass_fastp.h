//
//  bypass_fastp.h
//  
//
//  Created by Eric Strobel on 7/18/25.
//

#ifndef bypass_fastp_h
#define bypass_fastp_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

#include "../../../variant_maker/make_barcodes.h"

#include "../../../utils/io_management.h"

/* bypass_fastp: bypass fastp and perform simple read processing during testdata analysis. useful for systems in which fastp is not easily installed. */
void bypass_fastp(char * fq1, char * fq2, FILE ** ifp);

#endif /* bypass_fastp_h */

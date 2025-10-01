//
//  call_fastp_TPROBE.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef call_fastp_TPROBE_h
#define call_fastp_TPROBE_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

/* call_fastp_TPROBE: performs a system call to initiate fastp preprocessing, opens resulting files*/
int call_fastp_TPROBE(char * fq1, char * fq2, FILE ** ifp, fastp_params prms);

#endif /* call_fastp_TPROBE_h */

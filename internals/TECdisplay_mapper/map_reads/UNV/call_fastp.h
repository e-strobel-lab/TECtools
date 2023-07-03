//
//  call_fastp.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef call_fastp_h
#define call_fastp_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

/* call_fastp: performs a system call to initiate fastp preprocessing, opens resulting files*/
int call_fastp(names * nm, fastp_params prms);

#endif /* call_fastp_h */

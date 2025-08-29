//
//  fastp_params.h
//  
//
//  Created by Eric Strobel on 6/27/25.
//

#ifndef fastp_params_h
#define fastp_params_h

#include <stdio.h>

/* fastp_params: structure containing parameters for fastp processing */
typedef struct fastp_params {
    char path[MAX_LINE];  //path to fastp executable
    int mode;             //processing mode (multi-length or single-length)
    int limit;            //limit on number of reads to process
} fastp_params;           //NOTE: processing mode is also used by functions unrelated to fastp

#endif /* fastp_params_h */

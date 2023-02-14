//
//  appnd_att.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef appnd_att_h
#define appnd_att_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"

/* appnd_att: append attribute code to end of read id in line 1 of fastq file */
// Note: this function is not currently used, but is being kept in case of future utility
int appnd_att(char * read1, char * read2, char * att);

#endif /* appnd_att_h */

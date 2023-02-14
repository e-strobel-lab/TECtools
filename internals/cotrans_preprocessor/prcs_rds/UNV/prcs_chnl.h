//
//  prcs_chnl.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#ifndef prcs_chnl_h
#define prcs_chnl_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"
#include "../../cotrans_preprocessor_structs.h"

/* prcs_chnl: identify channel code from read 1 ID.
 these barcodes are used to split reads into modified and untreated channels
 barcodes:
 modified  = RRRYY
 untreated = YYYRR
 */
int prcs_chnl(char * read1_ID, metrics  * met);

#endif /* prcs_chnl_h */


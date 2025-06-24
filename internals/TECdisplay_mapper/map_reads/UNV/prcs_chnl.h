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

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

/* prcs_chnl: identify channel code from read 1 ID.
 these barcodes are used to split reads into modified and untreated channels
 barcodes:
 modified  = RRRYY
 untreated = YYYRR
 */
int prcs_chnl(char * read1_ID, TDSPLY_metrics  * met, int * chnl_match_type);

#endif /* prcs_chnl_h */


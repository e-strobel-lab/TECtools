//
//  filter_variants.h
//  
//
//  Created by Eric Strobel on 5/10/23.
//

#ifndef filter_variants_h
#define filter_variants_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/ispair.h"
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

/* filter_variants: filter variants that do not meet bases pair constraints*/
int filter_variants(char * variant, basemap * bmap);

#endif /* filter_variants_h */

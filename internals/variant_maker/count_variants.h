//
//  count_variants.h
//  
//
//  Created by Eric Strobel on 8/4/22.
//

#ifndef count_variants_h
#define count_variants_h

#include <stdio.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/gen_utils.h"
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

/* count_variants: determine the number of variants specified by the supplied variant templates */
int count_variants(basemap * bmap, int tot_templates, int mode);

#endif /* count_variants_h */

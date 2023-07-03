//
//  expand_variant_template.h
//  
//
//  Created by Eric Strobel on 8/3/22.
//

#ifndef expand_variant_template_h
#define expand_variant_template_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#include "../seq_utils/isIUPACbase.h"

#include "./mk_variant_nm.h"
#include "./make_barcodes.h"

/* expand_variant_template: recursively construct targets from a variant template base map */
void expand_variant_template(basemap * bmap, int nxt, char outpt[MAXLEN+1], struct trgt_vb lcl_bases, int mode);

#endif /* expand_variant_template_h */

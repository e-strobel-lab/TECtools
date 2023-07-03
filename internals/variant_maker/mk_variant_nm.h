//
//  mk_variant_nm.h
//  
//
//  Created by Eric Strobel on 8/3/22.
//

#ifndef mk_variant_nm_h
#define mk_variant_nm_h

#include <stdio.h>
#include <string.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#include "../seq_utils/basemap.h"

/* mk_variant_name: construct variant name from a local bases array */
void mk_variant_name(char * vrnt_nm, struct trgt_vb lcl_bases, struct basemap * bmap);

#endif /* mk_variant_nm_h */

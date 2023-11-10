//
//  parse_reference.h
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#ifndef parse_reference_h
#define parse_reference_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"
#include "../seq_utils/isIUPACbase.h"
#include "../seq_utils/isDNAbase.h"
#include "../seq_utils/basemap.h"

#include "../variant_maker/variant_maker_defs.h"
#include "../variant_maker/variant_maker_structs.h"
#include "../variant_maker/read_varFile.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

/* parse_reference: read reference lines in constraints file and store the reference sequence, variable base
 reference sequence, and parsed vbases in wt_source and basemap structs. parsing vbases from the variable base
 reference sequence establishes the reference that vbases in each constraint and variant names will be checked
 against to establish their validity. */
int parse_reference(FILE * ipt, basemap * bmap, wt_source * wt, char ** cnstnt_indels, int print);

#endif /* parse_reference_h */

//
//  search_4_vbase_match.h
//  
//
//  Created by Eric Strobel on 6/8/23.
//

#ifndef search_4_vbase_match_h
#define search_4_vbase_match_h

#include <stdio.h>

#include "../global/global_defs.h"

#include "../seq_utils/is_dgnrt_mtch.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

#define EXACT_VBASE_MATCH 0
#define DGNRT_VBASE_MATCH 1

/* search_for_vbase_match: search vbase haystack for match to vbase needle. if match
 is found, return pointer to needle match in haystack, otherwise return NULL */
base_params * search_4_vbase_match(base_params * haystack, base_params * needle, int vb_cnt, int mode);

#endif /* search_4_vbase_match_h */

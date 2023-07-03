//
//  search_4_vbase_match.c
//  
//
//  Created by Eric Strobel on 6/8/23.
//

#include <stdio.h>

#include "../global/global_defs.h"

#include "../seq_utils/is_dgnrt_mtch.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

#include "search_4_vbase_match.h"

/* search_for_vbase_match: search vbase haystack for match to vbase needle. if match
 is found, return pointer to needle match in haystack, otherwise return NULL */
base_params * search_4_vbase_match(base_params * haystack, base_params * needle, int vb_cnt, int mode)
{
    int i = 0; //general purpose index
    
    for (i = 0; i < vb_cnt; i++) {            //for every variable base in haystack
        if (needle->typ == haystack[i].typ && //check type match
            needle->pos == haystack[i].pos && //check position match
            needle->ins == haystack[i].ins) { //check insertion contig match
            
            switch (mode) {
                    
                case EXACT_VBASE_MATCH: //test for exact sequence match
                    if (needle->seq == haystack[i].seq) {
                        return &haystack[i];
                    }
                    break;
                 
                case DGNRT_VBASE_MATCH: //test for degenerate sequence match
                    if (is_dgnrt_mtch(needle->seq, haystack[i].seq)) {
                        return &haystack[i];
                    }
                    break;
                    
                default:
                    printf("search_for_vbase_match: error - unrecognized mode. aborting...");
                    abort();
                    break;
            }
        }
    }
    
    return NULL; //if no matches in haystack, return null pointer
}

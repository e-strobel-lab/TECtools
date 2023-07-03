//
//  filter_variants.c
//  
//
//  Created by Eric Strobel on 5/10/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../seq_utils/ispair.h"
#include "../seq_utils/basemap.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#include "filter_variants.h"

int allow_GU_base_pairs = 0; //flag to allow GU base pairs

/* filter_variants: filter variants that do not meet bases pair constraints*/
int filter_variants(char * variant, basemap * bmap)
{
    extern int allow_GU_base_pairs; //flag to allow GU base pairs
    
    int p = 0; //pair constraint index
    int i = 0; //general purpose index
    
    int min_pair = 0; //minimum base pair identity to allow
    
    if (allow_GU_base_pairs) { //if flag to allow GU base pairs is on
        min_pair = GU_PAIR;    //set minimum pair to GU
    } else {
        min_pair = AT_PAIR;    //otherwise set minimum pair to AT
    }
    
    for (p = 0; bmap->rP[p][0]; p++) {                                         //for every secondary structure constraint
        for (i = 0; variant[i] && i < MAXLEN; i++) {
            if (bmap->rP[p][i] == '(') {                                       //when the first mate of a base pair is found
                if (ispair(variant[i], variant[bmap->prs[p][i]]) < min_pair) { //check if it can pair with the second mate
                    return 1;                                                  //if not, return 1
                }
            }
        }
    }

    return 0; //return 0 if the variant passes filter
}

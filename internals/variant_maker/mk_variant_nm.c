//
//  mk_variant_nm.c
//  
//
//  Created by Eric Strobel on 8/3/22.
//

#include <stdio.h>
#include <string.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./variant_maker_defs.h"
#include "./variant_maker_structs.h"

#include "../seq_utils/basemap.h"

#include "mk_variant_nm.h"

/* mk_variant_name: construct variant name from a local bases array */
void mk_variant_name(char * vrnt_nm, struct trgt_vb lcl_bases, struct basemap * bmap)
{
    int i = 0;  //general purpose index
    int d = 0;  //deletion index
    
    char tmp[MAX_VBASE_NAME_FIELD+1]={0}; //array for assembling variant attributes during name construction
    
    for (i = 0; i < lcl_bases.cnt; i++) { //for every variable base that was identified
                
        mk_vbase_nm(bmap, lcl_bases.ix[i], tmp, MAX_VBASE_NAME_FIELD+1, lcl_bases.bs[i]); //assemble vbase name
        strcat(vrnt_nm, tmp); //append variable base to variant name
        
        if (i+1 < lcl_bases.cnt) {  //if there is another variable base
            strcat(vrnt_nm, "_");   //print an underscore separator
        }
    }

    return;
}

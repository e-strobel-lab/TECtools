//
//  parse_id_link.h
//  
//
//  Created by Eric Strobel on 9/25/26.
//

#ifndef parse_id_link_h
#define parse_id_link_h

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "../../TECdisplay_mapper/TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper/TECdisplay_mapper_structs.h"

#include "../../utils/io_management.h"

#define MAX_ID_LINK_FIELDS 4

/* var_id_link: storage for linked concise/verbose ids and associated barcodes */
typedef struct var_id_link {
    char * v;
    char * c;
    char * bc1;
    char * bc2;
    int bc_cnt;
    int subbed;
} var_id_link;

/* parse_id_link: parse id link file and store linked concise/verbose ids and barcode sequences in var_id_link structure */
void parse_id_link(var_id_link * id_link, int cnt, char * f_nm);

#endif /* parse_id_link_h */

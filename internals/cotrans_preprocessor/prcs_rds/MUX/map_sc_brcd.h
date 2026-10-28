//
//  map_sc_brcd.h
//  
//
//  Created by Eric Strobel on 10/23/26.
//

#ifndef map_sc_brcd_h
#define map_sc_brcd_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../seq_utils/C3SCbc_tree.h"

#include "../../../variant_maker/make_barcodes.h"

/* map_sc_brcd: map structure cassette barcode using binary tree search */
C3SC_barcode_node * map_sc_brcd(char * rd_ipt, C3SC_barcode_node * scbt_root);

#endif /* map_sc_brcd_h */

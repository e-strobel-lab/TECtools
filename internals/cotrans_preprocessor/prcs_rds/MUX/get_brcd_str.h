//
//  get_brcd_str.h
//  
//
//  Created by Eric Strobel on 10/17/25.
//

#ifndef get_brcd_str_h
#define get_brcd_str_h

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../variant_maker/make_barcodes.h"

/* get_brcd_str: get barcode string from UMI in read ID */
void get_brcd_str(char * brcd_str, char * read1_ID);

#endif /* get_brcd_str_h */

//
//  get_brcd_str.c
//  
//
//  Created by Eric Strobel on 10/17/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../variant_maker/make_barcodes.h"

#include "get_brcd_str.h"


/* get_brcd_str: get barcode string from UMI in read ID */
void get_brcd_str(char * brcd_str, char * read1_ID)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    //iterate to first space in read 1 name, which indicates end of channel barcode
    for (i = 0; read1_ID[i] != ' ' && read1_ID[i]; i++) { ;}
    
    //check that read1_ID contained space and that expected UMI location is within array bounds
    if (!read1_ID[i]) { //check that loop did not exit on null character
        printf("prcs_chnl: error - unexpected read1 id line format. aborting...\n");
        abort();
    } else if (i <= ((MAX_BARCODE_LEN*2)+1)) { //check that negative index won't go outside of array bounds
        printf("prcs_chnl: error - unexpected short read1 id line. aborting...\n");
        abort();
    }
    
    //set index to start of barcode string, which is 2 barcode lengths + 1 (to account for an '_') upstream of the space
    i -= ((MAX_BARCODE_LEN*2)+1);
    
    if (read1_ID[i-1] != ':') { //if preceding character is not a colon, throw error and abort
        printf("get_brcd_str: error - unexpected format for trimmed read 1/2 head sequences. aborting...\n");
        abort();
    }
    
    //store barcode string
    for (j = 0; j < MAX_BARCODE_LEN && read1_ID[i] != '_'; i++, j++) {
        brcd_str[j] = read1_ID[i];
    }
    brcd_str[j] = '\0';
    
    if (read1_ID[i] != '_' || j != MAX_BARCODE_LEN) { //check that barcode was copied as expected
        printf("get_brcd_str: error - unexpected format for trimmed read 1/2 head sequences. aborting...\n");
        abort();
    }
    
    return;
}

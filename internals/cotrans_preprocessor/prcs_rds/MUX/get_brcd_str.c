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
void get_brcd_str(char * brcd_rd1, char * brcd_rd2, int get_rd1_bc, int get_rd2_bc, char * read1_ID)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int spc_ix = 0; //index of the space follows the barcode sequences
        
    //iterate to first space in read 1 name, which indicates end of the barcode
    for (spc_ix = 0; read1_ID[spc_ix] != ' ' && read1_ID[spc_ix]; spc_ix++) { ;}
    
    //check that read1_ID contained space and that expected UMI location is within array bounds
    if (!read1_ID[spc_ix]) { //check that loop did not exit on null character
        printf("prcs_chnl: error - unexpected read1 id line format. aborting...\n");
        abort();
    } else if (spc_ix <= ((MAX_BARCODE_LEN*2)+1)) { //check that negative index won't go outside of array bounds
        printf("prcs_chnl: error - unexpected short read1 id line. aborting...\n");
        abort();
    }
    
    //get read 1 barcode
    if (get_rd1_bc) {
        
        //set index to start of bc string, which is 2 bc lengths + 1 (to account for an '_') upstream of the space
        i = spc_ix - ((MAX_BARCODE_LEN*2)+1);
        
        if (read1_ID[i-1] != ':') { //if preceding character is not a colon, throw error and abort
            printf("get_brcd_str: error - unexpected format for trimmed read 1/2 head sequences. aborting...\n");
            abort();
        }
        
        //store barcode string
        for (j = 0; j < MAX_BARCODE_LEN && read1_ID[i] != '_' && read1_ID[i]; i++, j++) {
            brcd_rd1[j] = read1_ID[i];
        }
        brcd_rd1[j] = '\0';
        
        //check loop end condition
        if (read1_ID[i] != '_' || j != MAX_BARCODE_LEN) { //check that barcode was copied as expected
            printf("get_brcd_str: error - unexpected format for trimmed read 1/2 head sequences. aborting...\n");
            abort();
        }
    }
    
    //get read 2 barcode
    if (get_rd2_bc) {
        
        //set index to start of bc string, which is 1 bc length upstream of the space
        i = spc_ix - MAX_BARCODE_LEN;
        
        if (read1_ID[i-1] != '_') { //if preceding character is not an underscore, throw error and abort
            printf("get_brcd_str: error - unexpected format for trimmed read 1/2 head sequences. aborting...\n");
            abort();
        }
        
        //store barcode string
        for (j = 0; j < MAX_BARCODE_LEN && read1_ID[i] != ' ' && read1_ID[i]; i++, j++) {
            brcd_rd2[j] = read1_ID[i];
        }
        brcd_rd2[j] = '\0';
        
        //check loop end condition
        if (read1_ID[i] != ' ' || j != MAX_BARCODE_LEN) { //check that barcode was copied as expected
            printf("get_brcd_str: error - unexpected format for trimmed read 1/2 head sequences. aborting...\n");
            abort();
        }
    }
    
    return;
}

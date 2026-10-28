//
//  map_sc_brcd.c
//  
//
//  Created by Eric Strobel on 10/23/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../cotrans_preprocessor_defs.h"

#include "../../../seq_utils/C3SCbc_tree.h"

#include "../../../variant_maker/constant_seqs.h"
#include "../../../variant_maker/make_barcodes.h"

#include "map_sc_brcd.h"

/* map_sc_brcd: map structure cassette barcode using binary tree search */
C3SC_barcode_node * map_sc_brcd(char * rd_ipt, C3SC_barcode_node * scbt_root)
{
    extern char c3_upper[20]; //uppercase C3 sequence
    extern char ttcg[5];      //TTCG tetraloop sequence
    
    //c3_rmdr is a pointer to the segment of C3 that remains after
    //processing by fastp, which trims 16 nts at the head of each read
    char * c3_rmdr = &c3_upper[MAX_BARCODE_LEN-VL_CHNL_BC_LEN];
    
    //allocate memory to store copy of read, then copy the read
    char * rd_cpy = NULL; //copy of input read, used for parsing
    if ((rd_cpy = malloc((strlen(rd_ipt)+1) * sizeof(*rd_cpy))) == NULL) {
        printf("C3SC_barcode_tree: error - failed to allocate memory for read 1 sequence copy during SC barcode mapping. aborting...\n");
        abort();
    }
    strcpy(rd_cpy, rd_ipt);
    
    char * p_c3_rmdr = NULL; //pointer to the C3 remainder sequence within the parsed target
    char * p_upstrm = NULL;  //pointer to the upstream segment of the SC barcode
    char * p_dnstrm = NULL;  //pointer to the downstream segment of the SC barcode
        
    C3SC_barcode_node * scbc = NULL; //pointer to mapped barcode node
        
    if ((p_c3_rmdr = strstr(rd_cpy, c3_rmdr)) != NULL) {                   //if a match to C3 remainder is found
        if ((uint64_t)(p_c3_rmdr) == (uint64_t)(rd_cpy)) {                 //and the match is at the head of the read
            p_upstrm = &p_c3_rmdr[strlen(c3_rmdr)];                        //set pointer to upstream SC barcode seg
            p_upstrm[MAX_C3SCBC_TREE_LAYER] = '\0';                        //terminate upstream barcode sequence
            scbc = srch_C3SCbc_tree(p_upstrm, 0, &p_upstrm[0], scbt_root); //determine if upstream BC maps to SC BC
            
            p_dnstrm = &p_upstrm[MAX_C3SCBC_TREE_LAYER + strlen(ttcg)];    //set ptr to calculated downstream BC start
            //TODO: need to add way to assess downstream barcode, however this will require data that informs how to search for the downstream barcode. The upstream barcode is straightforward to find because a primer masks any mutations in the C3 sequence in cDNA so that, in most cases, the start of the usptream barcode will immediately follow the C3 sequence
        }
    }
    
    free(rd_cpy);  //free allocated memory
    rd_cpy = NULL; //reset pointer to null
    
    return scbc; //return pointer to SC barcode node
}

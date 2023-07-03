//
//  variant_maker_structs.h
//  
//
//  Created by Eric Strobel on 8/3/22.
//

#ifndef variant_maker_structs_h
#define variant_maker_structs_h

#include <stdio.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./variant_maker_defs.h"

#include "../seq_utils/seq2bin_hash.h"

typedef struct fasta {   //structure to store fasta entries
    char * nm;           //sequence name
    char * sq;           //sequence with multiple sequence alignment characters preserved
} fasta;

typedef struct names {   //structure to store input names
    char iptV[MAX_LINE+1]; //input variant template file name
    char iptB[MAX_LINE+1]; //input barcode sequence file name
    char vTmp[MAX_LINE+1]; //variant template name
    char brcd[MAX_LINE+1]; //barcodes name
} names;

struct trgt_vb { //track variable base identity through recursion paths
    int cnt;                           //current index (== variable bases previously observed)
    char bs[SEQ2BIN_MAX_KEY+1];        //variable base identity
    int  ix[SEQ2BIN_MAX_KEY+1];        //variable base index, used to look up information in basemap struct
};


#endif /* variant_maker_structs_h */

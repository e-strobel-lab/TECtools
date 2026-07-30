//
//  get_msa_subseq.c
//  
//
//  Created by Eric Strobel on 11/24/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../seq_utils/isIUPACbase.h"

#include "get_msa_subseq.h"

/* get_msa_subseq: get subsequence from multiple sequence alignment line */
int get_msa_subseq(char ** seq, char * msa, int b1, int b2, int mode)
{
    //NOTE: for nucleotide identity descriptors, bound 1 is the index at
    //which the string starts and bound 2 is the length of the string.
    //for all other descriptors, bounds 1 and 2 are both indices.
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    char tmp_seq[MAX_LINE+1] = {0}; //array to temporarily store msa subsequence
    
    int idx_lmt = 0; //index limit
    int len_lmt = 0; //length limit
    
    //set limits based on mode. this determines whether index value or
    //sequence length will be used to limit the loop below
    if (mode == BOUND2_INDEX) { //if bound 2 is an index
        idx_lmt = b2;           //set idx_limit to bound 2
        len_lmt = MAX_LINE;     //set len_lmt to MAX_LINE length
        
    } else if (mode == BOUND2_LENGTH) { //if bound 2 is a length
        idx_lmt = MAX_LINE;             //set idx_limit to MAX_LINE length
        len_lmt = b2;                   //set len_limit to bound 2
    }
        
    //iterate from bound 1 index to the limit set by the mode and store the subsequence
    for (i = b1, j = 0; i <= idx_lmt && j < len_lmt && msa[i] && i < MAX_LINE; i++) {
        
        if (isIUPACbase(msa[i])) {                //only store IUPAC DNA/RNA seq chars, not spacers
            if (msa[i] == 't' || msa[i] == 'T') { //if t/T base,
                tmp_seq[j++] = 'U';               //convert to U
            } else {
                tmp_seq[j++] = toupper(msa[i]);   //otherwise, store base as uppercase
            }
        }
    }
    tmp_seq[j] = '\0'; //append terminating null char
    
    //allocate memory for storing the subsequence
    if ((*seq = malloc((j+1) * sizeof(**seq))) == NULL) {
        printf("get_msa_subseq: error - failed to allocate memory for subsequence storage. aborting...\n");
        abort();
    }
    
    strcpy(*seq, tmp_seq); //store msa subsequence
        
    return j; //return subsequence length
}

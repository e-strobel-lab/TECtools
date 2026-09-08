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
//NOTE: input sequence indexing must match the true numbering of the sequence (index 0 should have a placeholder char)
int get_msa_subseq(char ** seq, char * msa, int b1, int b2, int mode, char * pre, char * nxt, int set_pre, int set_nxt)
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
    
    } else {
        printf("get_msa_subseq: error unexpected mode\n");
        abort();
    }
    
    int fnd_bs = 0; //flag that a base was found when searching for prev/next base relative to the subseq window
    
    //search for and set previous/next bases
    if (set_pre || set_nxt) {
        
        if (mode == BOUND2_INDEX) { //pre/next base identification is only available for index bounds currently
            
            if (set_pre) { //set preceeding nucleotide
                for (i = b1-1, fnd_bs = 0; i >=1 && !fnd_bs; i--) { //search upstream for base character
                    if (isIUPACbase(msa[i])) {                      //if an IUPAC base is found
                        *pre = ret_upper_rna_nt(msa[i]);            //set pre as uppercase RNA nt
                        fnd_bs = 1;                                 //set fnd_bs flag to true
                    }
                }
                if (!fnd_bs) {   //if no base was found before the start of the string was reached
                    *pre = '\0'; //set pre to 0
                }
            }
            
            
            if (set_nxt) { //set next nucleotide
                for (i = b2+1, fnd_bs = 0; msa[i] && !fnd_bs; i++) { //search downstream for base character
                    if (isIUPACbase(msa[i])) {                       //if an IUPAC base is found
                        *nxt = ret_upper_rna_nt(msa[i]);             //set nxt as uppercase RNA nt
                        fnd_bs = 1;                                  //set fnd_bs flag to true
                    }
                }
                if (!fnd_bs) {   //if no base was found before the end of the string was reached
                    *nxt = '\0'; //set nxt to 0
                }
            }
        } else if (mode == BOUND2_LENGTH) {
            printf("get_msa_subseq: set_flanking is not currently enabled for BOUND2_LENGTH mode. aborting...\n");
            abort();
        }
    }
        
    //iterate from bound 1 index to the limit set by the mode and store the subsequence
    //skip all non-IUPAC base characters
    for (i = b1, j = 0; i <= idx_lmt && j < len_lmt && msa[i] && i < MAX_LINE; i++) {
        if (isIUPACbase(msa[i])) {
            tmp_seq[j++] = ret_upper_rna_nt(msa[i]); //set uppercase RNA nt
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

/* ret_upper_rna_nt: return uppercase RNA nucleotide base of input char */
char ret_upper_rna_nt(char c)
{
    if (c == 't' || c == 'T') {  //if t/T base,
        return 'U';              //convert to U
    } else if (isIUPACbase(c)) {
        return toupper(c);       //otherwise, store base as uppercase
    } else {
        printf("ret_upper_rna_nt: error - %c is not an RNA nucleotide. aborting...\n", c);
        abort();
    }
}

//
//  seq2bin_long.c
//  
//
//  Created by Eric Strobel on 3/6/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "./seq2bin_hash.h"
#include "seq2bin_long.h"

/* seq2bin_long: convert DNA sequence to two bit notation */
int seq2bin_long(char * hash_seq, binary_seq * bin_seq, int array_max) {
    
    extern int debug_S2B_hash;    //flag for running debug mode
    extern int trace_hash;        //flag to print 2-bit DNA seq encoding process
    
    debug_S2B_hash = 0;
    trace_hash = 0;
    
    //check that bin_seq pointer has been set
    if (bin_seq == NULL) {
        printf("seq2bin_long: error - bin_seq pointer is null.aborting...\n");
        abort();
    }
    
    int nts_per_block = (sizeof(bin_seq->sq[0])*8)/2; //number of nucleotides per bin_seq->sq index
    int filled_blocks = 0;  //number of blocks that will be completely filled by input sequence
    int partial_blocks = 0; //number of blocks that will be partially filled by input sequence

    bin_seq->ln = strlen(hash_seq); //set sequence length
    
    //check that 2-bit encoded sequences fits within the maximum-sized array
    if ((array_max * nts_per_block) < bin_seq->ln) {
        printf("seq2bin_long: error - insufficient array size (%d) for two-bit encoded sequence with length %lu nt. aborting...\n", array_max, strlen(hash_seq));
        abort();
    } else {
        filled_blocks = bin_seq->ln/nts_per_block;              //determine number of completely filled blocks
        partial_blocks = (bin_seq->ln % nts_per_block) ? 1 : 0; //determine number of partially filled blocks
        bin_seq->mx = filled_blocks + partial_blocks;           //max index is filled + partial blocks
        
        //allocate memory for storing 2-bit encoded sequence
        if ((bin_seq->sq = calloc(bin_seq->mx, sizeof(*bin_seq->sq))) == NULL) {
            printf("seq2bin_long: error - sequence memory allocation failed\n");
            return 1;
        }
    }

    int block = 0;     //block index
    int i = 0;         //general purpose index
    int mshft = 0;     //maximum bit shift value
    int shft = 0;      //bit shift value for each iteration of loop
    uint64_t ctob = 0; //value of nucleotide character to 2 bit conversion
    
    if (trace_hash) { printf("stepwise 2-bit DNA sequence encoding:\n");}
    
    //perform hash. in each iteration of the for loop a DNA sequence character is
    //converted to 2bit notation by masking the lowest 3 bits and rightshifting 1 bit.
    //this conversion yields a unique value for each DNA base (ACGT) as shown below.
    //
    //          vv-bits used for hash
    //A 65 01000001 00
    //C 67 01000011 01
    //G 71 01000111 11
    //T 84 01010100 10
    //              ^^-two bit encoding
    //
    //the 2bit-encoded DNA base is then leftshifted by 'shft' positions and OR-ed with
    //the 'bin_seq[block]' variable, which stores the 2-bit encoded DNA sequence in a uint64_t.
    //'shft' is initialized to 'mshft', which is the maximum bitshift minus 2 (e.g. ctob
    //is a 64 bit unsigned int so mshft=62. this shifts the first 2-bit encoded DNA base
    //into the highest two bits in the first iteration of the for loop, and 'shft' is
    //decremented by two after each iteration so that the next 2-bit encoded DNA base will
    //be leftshifted by two fewer bits that the previous 2-bit encoded DNA base. iterating
    //this process encodes the entire DNA sequence in 2bit notation within a uint64_t.
    //
    //it is possible for a DNA sequencing read to contain an 'N' base. N bases yield the
    //same 2bit value as G bases. if these reads yield a hash value for a node that has
    //been assigned a target, the read will be discarded at the string comparison step
    //in the srch_htbl function.
    
    mshft = (sizeof(bin_seq->sq[0]) << 3) - 2; //initialize maximum shift to number of bits - 2
    
    for (i = 0, block = 0; block < bin_seq->mx && hash_seq[i]; block++) {
        
        //test if input sequence character is a non-native DNA base
        if (hash_seq[i] != 'A' && hash_seq[i] != 'a' &&
            hash_seq[i] != 'T' && hash_seq[i] != 't' &&
            hash_seq[i] != 'G' && hash_seq[i] != 'g' &&
            hash_seq[i] != 'C' && hash_seq[i] != 'c') {
            bin_seq->nn = 1;
        }
        
        if (trace_hash) {                //if trace_hash is on...
            printf("block %d\n", block); //...print block index
        }
        
        
        for (ctob = 0, shft = mshft; hash_seq[i] && shft >= 0; i++, shft -= 2) {
            ctob = (hash_seq[i] & 7) >> 1;  //lowest 3 bit mask + 1 bit rshift produces unique int for ATGC
            bin_seq->sq[block] |= ctob << shft; //leftshift to OR 2bit encoded nucleotide into place
            
            //if trace_hash is on, print stepwise visualization
            //of 2bit DNA sequence encoding (not typically used)
            if (trace_hash) {
                printf("%c->%llu\t", hash_seq[i], (long long unsigned int)ctob);
                printbin(bin_seq->sq[block]);
            }
        }
        
        if (trace_hash) {
            printf("\n\n");
        }
    }
    
    //check that expected number of blocks was used when converting input sequence to 2-bit notation
    if (block != bin_seq->mx) {
        printf("seq2bin_long: error - number of blocks used to encode sequence does not match calculated value");
    }
    
    //in debug mode, print alignment of DNAseq -> 2bit notation conversion
    int j = 0;
    if (debug_S2B_hash) {
        for (i = 0, block = 0; block < bin_seq->mx && hash_seq[i]; block++) {
            printf("char%d:\t", block);
            for (j = 0; j < nts_per_block && hash_seq[i]; i++, j++) {
                printf(" %c", hash_seq[i]);
                if (!((i+1)%4)) {
                    printf(" ");
                }
            }
            printf("\n2bit%d:\t", block);
            printbin(bin_seq->sq[block]);
            printf("\n");
        }
    }
    
    return 1;
}

/* srch_htbl_prehash: search hash table for sequence match */
struct compact_h_node** srch_ctrg_htbl(binary_seq * bin_seq, uint64_t hash, compact_h_node **htbl, int trace_search)
{
    if (trace_search) {
        printf("in ctrg htbl\n");
    }

    
    compact_h_node **p_nd = NULL;  //pointer for h_node pointer handling
    p_nd = &htbl[hash];            //initialize node pointer to hash table entry
    
    uint64_t bytes2cmpr = sizeof(bin_seq->sq[0]) * bin_seq->mx; //number of bytes to compare
        
    //search hash table for match to input sequence
    //performs string comparison to check for true match
    //if not true match, proceed to next node in linked list and test again
    
    while ((*p_nd) != NULL) {  //while node has been assigned
        
        //compare input with hash table entry
        //both a bit comparison and a length comparison are performed. The length
        //comparison is necessary because A nucleotides are encoded as '00', so a
        //sequence with trailing A's can appear shorter than it is
        
        if (!memcmp(bin_seq->sq, (*p_nd)->ctrg->bsq.sq, bytes2cmpr) && bin_seq->ln == (*p_nd)->ctrg->bsq.ln) {
            
            if (trace_search) {
                printf("match, leaving ctrg htbl\n");
            }
            
            return p_nd;          //found match, return pointer to node
        } else {
            p_nd = &(*p_nd)->nxt; //proceed to next node in linked list
        }
    }
    
    if (trace_search) {
        printf("no match, leaving ctrg htbl\n");
    }
    
    return p_nd; //no match, return address to last node that was checked
}

/* print_bin_seq: print 2-bit encoded sequence */
void print_bin_seq(char * ipt)
{
    int i = 0;
    binary_seq bin_seq = {0};
    
    seq2bin_long(ipt, &bin_seq, MAX_BLOCKS);
    
    for (i = 0; i < bin_seq.mx; i++) {
        printf("%llu", (long long unsigned int)(bin_seq.sq[i]));
        if (i+1 == bin_seq.mx) {
            printf("\n");
        } else {
            printf("\t");
        }
    }
    
}

/* fprint_bin_seq: print 2-bit encoded sequence to file */
void fprint_bin_seq(FILE * out_fp, char * ipt)
{
    int i = 0;
    binary_seq bin_seq = {0};
    
    seq2bin_long(ipt, &bin_seq, MAX_BLOCKS);
    
    for (i = 0; i < bin_seq.mx; i++) {
        fprintf(out_fp, "%llu", (long long unsigned int)(bin_seq.sq[i]));
        if (i+1 == bin_seq.mx) {
            fprintf(out_fp, "\n");
        } else {
            fprintf(out_fp, "\t");
        }
    }
}



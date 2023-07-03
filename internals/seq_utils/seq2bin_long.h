//
//  seq2bin_long.h
//  
//
//  Created by Eric Strobel on 3/6/23.
//

#ifndef seq2bin_long_h
#define seq2bin_long_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "./seq2bin_hash.h"

#define MAX_BLOCKS 256
#define NTS_PER_BLOCK 32
#define MAX_DEPTH_TRACKED 1024

/* binary_seq: structure for storing 2-bit encoded DNA sequences */
typedef struct binary_seq {
    uint64_t * sq; //2-bit encoded sequence
    uint8_t ln;    //number of nucleotides in 2-bit encoded sequence
    uint8_t mx;    //max index of 2-bit encoded sequence
    uint8_t nn;    //flag that sequence contains non-native DNA character
} binary_seq;

/* compact_target: memory efficient target structure*/
//TODO: COMPRESS MUL AND BL INTO ONE UINT8_T
typedef struct compact_target {
    uint64_t bid;    //target identifier //TODO: consider renaming to 'tid' for consistency
    binary_seq bsq;  //2-bit encoded target sequence
    uint8_t mul;     //flag that sequence is identical to a prior target
    uint8_t bl;      //flag that target is blacklisted
    uint32_t cnt;    //number of reads that map to target
    void * opt;
} compact_target;

/* compact_h_node: hash table node for compact targets */
typedef struct compact_h_node {
    compact_target *ctrg;        //pointer to target structure
    struct compact_h_node *nxt;  //points to next node in linked list, used to handle collisions
} compact_h_node;

/* compact_h_node_bank: hash table node bank for memory management */
typedef struct compact_h_node_bank {
    struct compact_h_node *chn;       //pointer for compact_h_node memory allocation
    struct compact_h_node_bank *nxt;  //pointer to next bank
    int count;                        //number of compact_h_nodes used in current bank
} compact_h_node_bank;


int seq2bin_long(char * hash_seq, binary_seq * bin_seq, int array_max);
struct compact_h_node** srch_ctrg_htbl(binary_seq * bin_seq, uint64_t hash, compact_h_node **htbl, int trace_search);
void print_bin_seq(char * ipt);
void fprint_bin_seq(FILE * out_fp, char * ipt);


#endif /* seq2bin_long_h */

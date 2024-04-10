//
//  seq2bin_hash.h
//  
//
//  Created by Eric Strobel on 3/15/22.
//
//  230207: memory for target struct members id, sq, and rc is now dynamically allocated in parse_3pEnd_trgts.c

#ifndef seq2bin_hash_h
#define seq2bin_hash_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define SEQ2BIN_MAX_IPT 2048    //maximum target length
#define SEQ2BIN_MAX_KEY 32		//maximum key length
#define TABLE_SIZE 16000057		//hash table size
#define BLOCK_SIZE 131072		//memory allocation block size

/* target: structure containing values for hash table targets */
typedef struct target {
    char *key;  //pointer to key for generating hash, can point to sq, rc, or opt member
    char * id;	//target identifier
    char * sq;	//target sequence
    char * rc;	//reverse complement of target sequence
    int cnt;    //number of reads that map to target
    int mul;    //flag that sequence is identical to a prior target
    void * opt;	//pointer for adding additional data to target structure
} target;

/*h_node: hash table node */
typedef struct h_node {
    target *trg;				//pointer to target structure
    struct h_node *nxt;			//points to next node in linked list, used to handle collisions
} h_node;

/* h_node_bank: hash table node bank for memory management */
typedef struct h_node_bank {
    struct h_node *hn;			//pointer for h_node memory allocation
    struct h_node_bank *nxt;	//pointer to next bank
    int count;					//number of h_nodes used in current bank
} h_node_bank;

/* srch_htbl: hash table search function. searches hash table (htbl) for match to query.
 If a match is found, a pointer to the node that matches the query is returned. If a match
 is not found, NULL is returned. */
struct h_node** srch_htbl(char *query, h_node **htbl);

/* seq2bin_hash: convert DNA seq to 2 bit notation and return <2 bit seq value> % TABLE_SIZE. */
uint64_t seq2bin_hash(char *hash_seq);

/* extend_h_bank: increase hash table node bank size */
void extend_h_bank(h_node_bank *crrnt_hn_bank);

/* printbin: print binary representation of uint64_t value */
void printbin(uint64_t x);

#endif /* seq2bin_hash_h */

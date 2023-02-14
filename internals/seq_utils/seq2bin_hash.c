//
//  seq2bin_hash.c
//  
//
//  Created by Eric Strobel on 3/15/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "seq2bin_hash.h"

int debug_S2B_hash = 0; //flag for running debug mode
int trace_hash = 0;

//todo add check that local max does not exceed seq2bin hash max field

/* srch_htbl: search hash table for sequence match */
struct h_node** srch_htbl(char *query, h_node **htbl)
{
    //check that query length is less than maximum key length
    if (strlen(query) > SEQ2BIN_MAX_KEY) {
        printf("srch_htbl: error - query string exceeds (%lu chars) max key length (%d chars). aborting...\n", strlen(query), SEQ2BIN_MAX_KEY);
        abort();
    }
    
    h_node **p_nd = NULL;				//pointer for h_node pointer handling
    uint64_t hash = 0xFFFFFFFFFFFFFFFF;	//hash value
    
    //hash table search
    hash = seq2bin_hash(query);	//generate hash value
    p_nd = &htbl[hash];			//initialize node pointer to hash table entry
        
    //search hash table for match to input sequence
    //performs string comparison to check for true match
    //if not true match, proceed to next node in linked list and test again
    while ((*p_nd) != NULL) {						//while node has been assigned
        if (!strcmp(query, (*p_nd)->trg->key)) {	//compare input string with hash table entry
            return p_nd;							//found match, return pointer to node
        } else {
            p_nd = &(*p_nd)->nxt;					//proceed to next node in linked list
        }
    }
    return p_nd; //no match, return address to last node that was checked
}



/* seq2bin_hash: convert DNA sequence to two bit notation */
uint64_t seq2bin_hash(char *hash_seq) {
    
    extern int debug_S2B_hash;	//flag for running debug mode
    extern int trace_hash;		//flag to print 2-bit DNA seq encoding process
    
    int i = 0;
    int mshft = 0;               //maximum bit shift value
    int shft = 0;                //bit shift value for each iteration of loop
    uint64_t ctob = 0;           //value of nucleotide character to 2 bit conversion
    uint64_t hash_outpt = 0;     //2 bit encoded nucleotide character sequence
    
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
    //the 'hash_outpt' variable, which stores the 2-bit encoded DNA sequence in a uint64_t.
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
    
    mshft = (sizeof(hash_outpt) << 3) - 2; //initialize maximum shift to number of bits - 2
    for (i = 0, ctob = 0, shft = mshft; hash_seq[i] && shft >= 0; i++, shft -= 2) {
        ctob = (hash_seq[i] & 7) >> 1;	//lowest 3 bit mask + 1 bit rshift produces unique int for ATGC
        hash_outpt |= ctob << shft;		//leftshift to OR 2bit encoded nucleotide into place
        
        //if trace_hash is on, print stepwise visualization
        //of 2bit DNA sequence encoding (not typically used)
        if (trace_hash) {
            printf("%c->%llu\t", hash_seq[i], (long long unsigned int)ctob);
            printbin(hash_outpt);
        }
    }
    
    //in debug mode, print alignment of DNAseq -> 2bit notation conversion
    if (debug_S2B_hash) {
        printf("char:\t");
        for (i = 0; hash_seq[i]; i++) {
            printf(" %c", hash_seq[i]);
            if (!((i+1)%4)) {
                printf(" ");
            }
        }
        printf("\n2bit:\t");
        printbin(hash_outpt);
    }
    
    return hash_outpt % TABLE_SIZE; //return 2bit seq value % TABLE_SIZE
}


/* extend_h_bank: used to increase the size of the has table node bank. Not currently used. */
void extend_h_bank(h_node_bank *crrnt_hn_bank)
{
    if ((crrnt_hn_bank->nxt = calloc(1, sizeof(*(crrnt_hn_bank->nxt)))) == NULL) {
        printf("srch_h_tbl: error - hash table node bank memory allocation failed\n");
        abort();
    }
    crrnt_hn_bank = crrnt_hn_bank->nxt;
    if ((crrnt_hn_bank->hn = calloc(BLOCK_SIZE, sizeof(*(crrnt_hn_bank->hn)))) == NULL) {
        printf("srch_h_tbl: error - hash table node bank memory allocation failed\n");
        abort();
    }
}


/* printbin: print 64 bit int in binary notation */
void printbin(uint64_t x)
{
    long int i = 0;		//used to decrement from msb to zero when testing bits for 1 or 0 value
    long int j = 0;		//used to count printed bits so every byte can be delimited by a space
    long int msb = 0;	//power value of the most significant bit. msb=63 for uint64_t
    uint64_t pow2 = 0;	//used to test bit value and zero ON bits
    
    msb = ((sizeof(x) << 3) - 1); //msb = <size of x in bytes> * 8 - 1
    pow2 = pow(2, msb);    //set pow2 so that only the most significant bit is on
    
    //starting at the most significant bit, each iteration of the loop tests
    //whether one bit of the 64 bit int x is on. this is accomplished as follows:
    //
    //pow2 is initialized so that only the msb is on
    //each iteration of the loop tests whether x-pow2 >= 0
    //if x-pow2 >= 0, the bit was on, print 1, and subtract pow2 from x to zero the bit.
    //if x-pow2 < 0, the bit was off, print 0
    //at the end of each iteration, pow2 is rightshifted 1 bit so that
    //the value of the next lowest bit can be tested.
    
    for (i = msb, j = 0; i >= 0; i--, pow2 = pow2 >> 1) {
        if (((int64_t)(x - pow2)) >= 0) {	//bit i is on
            x -= pow2;		//subtract pow2 to zero bit i
            putchar('1');	//print 1
        } else {			//bit i is off
            putchar('0');	//print 0
        }
        if (j++ == 7) {		//print a space every 8 iterations
            putchar(' ');	//to separate bytes
            j = 0;			//reset counter
        }
    }
    putchar('\n');	//print newline
    
    return;
}

//
//  mk_dot_bracket_htbl.c
//  
//
//  Created by Eric Strobel on 6/30/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "../../seq_utils/seq2bin_hash.h"
#include "../../seq_utils/seq2bin_long.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"
#include "../set_seq_attributes/parse_descriptor_file.h"

#include "mk_dot_bracket_htbl.h"

/* mk_dot_bracket_htbl: makes compact target hash table */
/* hash table has linked list buckets for possible collisions */
void mk_dot_bracket_htbl(compact_h_node ** htbl, compact_h_node_bank * bank, compact_target * ctrg, int ctrg_cnt, structProps ** sp_list/*, target_params * trg_prms*/)
{
    compact_h_node **p_rdnd = NULL; //pointer for h_node handling
    
    int i = 0;              //general purpose index
    int j = 0;              //general purpose index
    int new_node = 0;       //counts number of nodes that are assigned a target
    int redundant = 0;      //counts number of redundant targets (same seq as prev target)
    int collisions = 0;     //counts number of collisions between barcodes
    
    int error = 0; //flag that error occurred and program should be aborted
    
    uint64_t hash = 0;  //hash key
    
    for (i = 0; i < ctrg_cnt; i++) {  //perform loop for each target
        
        hash = hash_long_bsq(&ctrg[i].bsq);                   //hash binary encoded sequence
        p_rdnd = srch_ctrg_htbl(&ctrg[i].bsq, hash, htbl, 0); //search hash table
        
        if ((*p_rdnd) == NULL) {                    //no existing hash table node for target sequence
            (*p_rdnd) = &bank->chn[bank->count++];  //assign node from hash node bank
            if (bank->count == BLOCK_SIZE) {        //check that bank was not filled
                extend_ch_bank(bank);               //extend bank if needed
            }
            (*p_rdnd)->ctrg = &(ctrg[i]);           //set node to point to current target
            (*p_rdnd)->ctrg->mul++;                 //track number of structProps that match first ctrg for current structure
            new_node++;                             //increment new_node counter
            
        } else {
            (*p_rdnd)->ctrg->mul++; //track number of structProps that match first ctrg for current structure
            ctrg[i].bl++;           //set flag that current target is redudant with previous target
            redundant++;            //increment redundant counter
        }
    }
    
    init_opt_db(ctrg, ctrg_cnt);                     //allocate opt_db memory in ctrg structures
    map_structProps_2_htbl(sp_list, ctrg_cnt, htbl); //
    
    //print hash table generation details
    printf("\n********************  hash table generation  ********************\n");
    printf("%6d dot-bracket structures were assessed\n", i);
    printf("%6d dot-bracket structures were assigned a node\n", new_node);
    printf("%6d dot-bracket structures were not assigned a node\n\n\n", redundant);
    
    //trg_prms->nr_cnt = new_node;  //set non-redundant node count
    //trg_prms->rdndnt = redundant; //set redundant node count
    
    opt_db * p_db_val = NULL;
    
    //print related structures
    for(i = 0; i < ctrg_cnt; i++) {
        if (!ctrg[i].bl) {
            p_db_val = ctrg[i].opt;
            printf("mul = %d, cnt = %d\n", ctrg[i].mul, ctrg[i].cnt);
            for (j = 0; j < ctrg[i].mul; j++) {
                printf("%03d: %s\n", j, p_db_val->sp[j]->db);
            }
            printf("\n\n");
        }
    }
}

/* dot_bracket_to_bin_hash: convert DNA sequence to two bit notation */
int dot_bracket_to_bin_hash(char * hash_seq, binary_seq * bin_seq, int array_max) {
    
    extern int debug_S2B_hash; //flag for running debug mode
    extern int trace_hash;     //flag to print 2-bit DNA seq encoding process
    
    debug_S2B_hash = 0;
    trace_hash = 0;
    
    //check that bin_seq pointer has been set
    if (bin_seq == NULL) {
        printf("dot_bracket_to_hash: error - bin_seq pointer is null. aborting...\n");
        abort();
    }
    
    int nts_per_block = (sizeof(bin_seq->sq[0])*8)/2; //number of nucleotides per bin_seq->sq index
    int filled_blocks = 0;  //number of blocks that will be completely filled by input sequence
    int partial_blocks = 0; //number of blocks that will be partially filled by input sequence

    bin_seq->ln = strlen(hash_seq); //set sequence length
    
    //check that 2-bit encoded sequences fits within the maximum-sized array
    if ((array_max * nts_per_block) < bin_seq->ln) {
        printf("dot_bracket_to_hash: error - insufficient array size (%d) for two-bit encoded sequence with length %lu nt. aborting...\n", array_max, strlen(hash_seq));
        abort();
    } else {
        filled_blocks = bin_seq->ln/nts_per_block;              //determine number of completely filled blocks
        partial_blocks = (bin_seq->ln % nts_per_block) ? 1 : 0; //determine number of partially filled blocks
        bin_seq->mx = filled_blocks + partial_blocks;           //max index is filled + partial blocks
        
        //allocate memory for storing 2-bit encoded sequence
        if ((bin_seq->sq = calloc(bin_seq->mx, sizeof(*bin_seq->sq))) == NULL) {
            printf("dot_bracket_to_hash: error - sequence memory allocation failed\n");
            abort();
        }
    }

    int block = 0;     //block index
    int i = 0;         //general purpose index
    int mshft = 0;     //maximum bit shift value
    int shft = 0;      //bit shift value for each iteration of loop
    uint64_t ctob = 0; //value of nucleotide character to 2 bit conversion
    
    if (trace_hash) { printf("stepwise 2-bit DNA sequence encoding:\n");}
    
    //perform hash. in each iteration of the for loop a dot-bracket sequence character
    //is converted to 2bit notation by masking all but the lowest 2 bits. this conversion
    //yields a unique value for each dot-bracket character as shown below.
    //
    //            vv-bits used for hash
    // ( 40 00101000 00
    // ) 41 00101001 01
    // . 46 00101110 10
    //               ^^-two bit encoding
    //
    //the 2bit-encoded character is then leftshifted by 'shft' positions and OR-ed with
    //the 'bin_seq[block]' variable, which stores the 2-bit encoded characters in a uint64_t.
    //'shft' is initialized to 'mshft', which is the maximum bitshift minus 2 (e.g. ctob
    //is a 64 bit unsigned int so mshft=62. this shifts the first 2-bit encoded DNA base
    //into the highest two bits in the first iteration of the for loop, and 'shft' is
    //decremented by two after each iteration so that the next 2-bit encoded DNA base will
    //be leftshifted by two fewer bits that the previous 2-bit encoded DNA base. iterating
    //this process encodes the entire DNA sequence in 2bit notation within a uint64_t.
    
    mshft = (sizeof(bin_seq->sq[0]) << 3) - 2; //initialize maximum shift to number of bits - 2
    
    for (i = 0, block = 0; block < bin_seq->mx && hash_seq[i]; block++) {
        
        if (trace_hash) {                //if trace_hash is on...
            printf("block %d\n", block); //...print block index
        }
        
        for (ctob = 0, shft = mshft; hash_seq[i] && shft >= 0; i++, shft -= 2) {
            
            //test if input sequence character is a non-dot-bracket character
            if (hash_seq[i] != '.' && hash_seq[i] != '(' && hash_seq[i] != ')') {
                printf("dot_bracket_to_hash: error - string contains unexpected character %c (ASCII %d). aborting...\n", hash_seq[i], hash_seq[i]);
            }
            
            ctob = (hash_seq[i] & 3);           //lowest 2 bit mask produces unique int for .() chars
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
        printf("dot_bracket_to_hash: error - number of blocks used to encode sequence does not match calculated value");
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

/* init_opt_db: allocate opt_db structure memory in compact_target structures */
void init_opt_db(compact_target * ctrg, int ctrg_cnt)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    opt_db * p_db_val = NULL; //pointer for dereferencing ctrg opt member as opt_db *
    
    for (i = 0; i < ctrg_cnt; i++) {             //for every target
        if (!ctrg[i].bl) {                       //if the target was not blacklisted (redundant with prev target)
            p_db_val = (opt_db *)(ctrg[i].opt);  //dereference opt as opt_db *
                                    
            //allocate memory for structProps structures. the mul member of ctrg, which is being used to count
            //the number of dot-bracket structures with the same secondary structure, is used to determine how
            //many structProps structures to allocate
            if ((p_db_val->sp = calloc(ctrg[i].mul, sizeof(*p_db_val->sp))) == NULL) {
                printf("init_opt_db: error - failed to allocate memory for structProps pointers. aborting...\n");
                abort();
            }
        }
    }
    
    return;
}

/* map_structProps_2_htbl: associate structProps structures with the relevant targets */
//TODO: should there be a test at the end that cnt == mul?
void map_structProps_2_htbl(structProps ** sp_list, int sp_cnt, compact_h_node ** htbl)
{
    int i = 0; //general purpose index
    
    compact_h_node ** p_rdnd = NULL; //pointer to hash table node pointer
    compact_h_node * chnNULL = NULL; //null compact_h_node pointer
    
    binary_seq bsq = {0};     //variable for storing binary-encoded barcode sequence
    uint64_t hash = 0;        //hash value
    opt_db * p_db_val = NULL; //pointer for dereferencing opt_db pointer
    
    for (i = 0; i < sp_cnt; i++) {                                 //for every structProps structure
        dot_bracket_to_bin_hash(sp_list[i]->db, &bsq, MAX_BLOCKS); //generate binary-encoded structure
        hash = hash_long_bsq(&bsq);                                //hash the binary-encoded structure
        p_rdnd = srch_ctrg_htbl(&bsq, hash, htbl, 0);              //search the hash table for a node that matches the structure
        
        if (*p_rdnd != NULL) {                                     //if a match was found
            p_db_val = (opt_db *)((*p_rdnd)->ctrg->opt);           //dereference opt as an opt_db_*
            p_db_val->sp[(*p_rdnd)->ctrg->cnt++] = sp_list[i];     //point the next structProps ptr to the current structProps
            
        } else {
            printf("map_structProps_2_htbl: srch_ctrg_htbl returned null node. aborting...\n");
            abort();
        }
    }
    
    return;
}

//
//  get_key.c
//  
//
//  Created by Eric Strobel on 9/30/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../../global/global_defs.h"
#include "../../../global/global_structs.h"

#include "../../../utils/debug.h"

#include "../../TECdisplay_mapper_defs.h"
#include "../../TECdisplay_mapper_structs.h"

#include "../../../seq_utils/seq2bin_hash.h"

#include "get_key.h"


/* get_key: generate key string composed the nucleotides at variable
 base positions in the input read sequence */
int get_key(char * key, char * end5p, char * qscore5p, char * minQv, target *refs, int key_type)
{
    extern int debug;
    
    int i = 0;                                     //general purpose index
    int len = strlen(end5p);                       //length of the input read sequence
    opt_ref * crnt_ref_val = (opt_ref *)refs->opt; //pointer to reference target optional values
    
    if (crnt_ref_val->vb_cnt > SEQ2BIN_MAX_KEY) {  //check that the number of variable bases is allowable
        printf("get_key: error - the number of variable bases specified by a reference sequence (%d) exceeds the maximum key length (%d). aborting...\n", crnt_ref_val->vb_cnt, SEQ2BIN_MAX_KEY);
    }
    
    //construct a hash table key using the variable base positions
    //specified by the reference sequence. checks are performed to
    //ensure that the read sequence is long enough to contain a
    //nucleotide at each variable base position, that the variable
    //bases are of sufficient quality, and that the resulting key
    //contains the number of bases specified by the reference sequence
    
    for (i = 0; i < crnt_ref_val->vb_cnt && i < SEQ2BIN_MAX_KEY; i++) {
 
        if (crnt_ref_val->vb_pos[i] < len) { // check that vb_pos index doesn't exceed string length
            
            if (key_type == TARGET_KEY) {                 //input sequence is a target, no qscore check needed
                key[i] = end5p[crnt_ref_val->vb_pos[i]];  //copy variable position base to key string
                
            } else if (key_type == READ_KEY && qscore5p[crnt_ref_val->vb_pos[i]] >= *minQv) {
                //input sequence is a read, check that qscore exceeds minimum qscore value
                key[i] = end5p[crnt_ref_val->vb_pos[i]];  //copy variable position base to key string
                
            } else {
                if (key_type != TARGET_KEY && key_type != READ_KEY) {
                    printf("get_key: error - unrecognized key type. aborting...\n");
                    abort();
                } else {
                    if (debug) {printf("variable base is too low quality to generate key\n");}
                    return 0; //minimum qscore not met, return 0
                }
            }
        } else {
            if (debug) {printf("read is too short to generate key\n");}
            return 0; //sequence too short to get all variable bases, return 0
        }
    }
    key[i] = '\0'; //terminate key string
    
    //check that key string contains the correct number of variable bases
    if ((i != crnt_ref_val->vb_cnt) || (strlen(key) != crnt_ref_val->vb_cnt)) {
        if (debug) {printf("key does not contain the correct number of variable bases\n");}
        return 0;
    } else {
        return 1;
    }
}

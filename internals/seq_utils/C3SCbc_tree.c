//
//  C3SCbc_tree.c
//  
//
//  Created by Eric Strobel on 10/20/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../variant_maker/constant_seqs.h"

#include "./revcomp.h"

#include "C3SCbc_tree.h"

char C3SCbc_node_left_key  = 'G'; //key for left C3SCbc node child
char C3SCbc_node_right_key = 'C'; //key for right C3SCbc node child

/* make_C3SCbc_tree: construct C3SC barcode tree */
void make_C3SCbc_tree(C3SC_barcode_node * scbc_tree, char key, int layer)
{
    int i = 0; //general purpose index
    
    static char crnt_str[MAX_C3SCBC_TREE_LAYER+1] = {0}; //current string
    
    if (layer == 0) { //if at layer zero, ensure all structure members are 0/NULL
        scbc_tree->key = '\0';
        scbc_tree->sq = NULL;
        scbc_tree->l  = NULL;
        scbc_tree->r  = NULL;
        scbc_tree->rc = NULL;
        
    } else { //otherwise, set key for current node and append key to current string
        scbc_tree->key = key;
        crnt_str[layer-1] = key;
    }
    
    
    if (layer < MAX_C3SCBC_TREE_LAYER) { //if the current layer is less than the max layer
        
        //allocate memory for child nodes
        if (((scbc_tree->l = calloc(1, sizeof(*scbc_tree->l))) == NULL) ||
            ((scbc_tree->r = calloc(1, sizeof(*scbc_tree->r))) == NULL)) {
            printf("make_C3SCbc_tree: error - failed to allocate memory for child nodes when generating tree. aborting...\n");
            abort();
        }
        
        //proceed to next layer for each child node
        make_C3SCbc_tree(scbc_tree->l, C3SCbc_node_left_key,  layer+1);
        make_C3SCbc_tree(scbc_tree->r, C3SCbc_node_right_key, layer+1);
        
    } else if (layer == MAX_C3SCBC_TREE_LAYER) { //if the current layer is the max layer
        
        //allocate memory for storing the sequence string
        if ((scbc_tree->sq = malloc((strlen(crnt_str)+1) * sizeof(*scbc_tree->sq))) == NULL) {
            printf("make_C3SCbc_tree: error - failed to allocate memory for sequence. aborting...\n");
            abort();
        }
        strcpy(scbc_tree->sq, crnt_str); //store the sequence string
        
    } else {
        printf("make_C3SCbc_tree: error - layer (%d) exceeded max layers (%d). aborting...\n", layer, MAX_C3SCBC_TREE_LAYER);
        abort();
    }
    
    if (layer == 0) { //if at layer 0 after generating tree, establish revcomp links
        set_C3SCbc_tree_revcomp(scbc_tree);
    }
    
    return;
}

/* srch_C3SCbc_tree: search for a match to str using C3SC barcode tree */
C3SC_barcode_node * srch_C3SCbc_tree(char * str, int i, char * nxt_key, C3SC_barcode_node * scbt_node)
{
    //navigate tree using keys from input string until the end of the input string is reached
    
    if (*nxt_key == C3SCbc_node_left_key) { //proceed to left child
        return (scbt_node = srch_C3SCbc_tree(str, i+1, &str[i+1], scbt_node->l));
        
    } else if (*nxt_key == C3SCbc_node_right_key) { //proceed to right child
        return (scbt_node = srch_C3SCbc_tree(str, i+1, &str[i+1], scbt_node->r));
        
    } else if (*nxt_key == '\0') { //reached end of input string
                 
        if (!strcmp(str, scbt_node->sq)) { //if search string matches node string,
            return scbt_node;              //return pointer to node
            
        } else { //failure should only happen if the node string is a substring of the search string beginning at index 0
            printf("srch_C3SCbc_tree: error - search string (%s) is not equivalent to matched node string (%s). is the node string a substring of the search string? aborting...\n", str, scbt_node->sq);
            abort();
        }
        
    } else if (*nxt_key == 'A' || *nxt_key == 'T' || *nxt_key == 'N'){
        return NULL; //A, T, and N may occur within a sequencing read but cause a failed search
        
    } else { //invalid character, throw error and abort
        printf("srch_C3SCbc_tree: invalid character '%c' in search string. aborting...\n", *nxt_key);
        abort();
    }
}

/* set_C3SCbc_tree_revcomp: link C3SC barcode nodes that are the reverse complement of each other */
void set_C3SCbc_tree_revcomp(C3SC_barcode_node * scbt_root)
{
    //list of possible barcodes
    char SCbrcds[C3SC_BRCD_CNT][MAX_C3SCBC_TREE_LAYER+1] = {"GGGG",
                                                            "GGGC",
                                                            "GGCG",
                                                            "GGCC",
                                                            "GCGG",
                                                            "GCGC",
                                                            "GCCG",
                                                            "GCCC",
                                                            "CGGG",
                                                            "CGGC",
                                                            "CGCG",
                                                            "CGCC",
                                                            "CCGG",
                                                            "CCGC",
                                                            "CCCG",
                                                            "CCCC"};

    
    //list of possible barcodes
    /*
    char SCbrcds[C3SC_BRCD_CNT][MAX_C3SCBC_TREE_LAYER+1] = {"gggg",
                                                            "gggc",
                                                            "ggcg",
                                                            "ggcc",
                                                            "gcgg",
                                                            "gcgc",
                                                            "gccg",
                                                            "gccc",
                                                            "cggg",
                                                            "cggc",
                                                            "cgcg",
                                                            "cgcc",
                                                            "ccgg",
                                                            "ccgc",
                                                            "cccg",
                                                            "cccc"};
    */
    
    int i = 0; //general purpose index
    
    char rc[MAX_C3SCBC_TREE_LAYER+1] = {0}; //reverse complement of barcode
    
    C3SC_barcode_node * sq_node = NULL; //current sequence node
    C3SC_barcode_node * rc_node = NULL; //reverse complement of current sequence node
    
    for (i = 0; i < C3SC_BRCD_CNT; i++) { //for every possible barcode
        
        reverse_complement(rc, SCbrcds[i], REVCOMP); //generate the reverse complement
        
        sq_node = srch_C3SCbc_tree(SCbrcds[i], 0, &SCbrcds[i][0], scbt_root); //identify the current sequence node
        rc_node = srch_C3SCbc_tree(rc, 0, &rc[0], scbt_root);                 //identify the revcomp node
        
        if (sq_node->rc == NULL && rc_node->rc == NULL) { //if both revcomp pointers are null, set the pointers
            sq_node->rc = rc_node;
            rc_node->rc = sq_node;
            
        } else if (sq_node->rc != NULL && rc_node->rc != NULL) { //do nothing, revcomp node pointers were already set
            ;
            
        } else { //throw error if set status of revcomp nodes is discordant
            printf("set_C3SCbc_tree_revcomp: error - discordant set status for revcomp pair. aborting...\n");
            abort();
        }
    }
}

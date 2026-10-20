//
//  C3SCbc_tree.h
//  
//
//  Created by Eric Strobel on 10/20/26.
//

#ifndef C3SCbc_tree_h
#define C3SCbc_tree_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../variant_maker/constant_seqs.h"

#include "./revcomp.h"

#define C3SCBC_KEY_CNT 2
#define MAX_C3SCBC_TREE_LAYER 4

/* C3SC_barcode_node: node for use in C3SC barcode tree */
typedef struct C3SC_barcode_node {
    char key;
    char * sq;
    struct C3SC_barcode_node * l;
    struct C3SC_barcode_node * r;
    struct C3SC_barcode_node * rc;
} C3SC_barcode_node;

/* make_C3SCbc_tree: construct C3SC barcode tree */
void make_C3SCbc_tree(C3SC_barcode_node * scbc_tree, char key, int layer);

/* srch_C3SCbc_tree: search for a match to str using C3SC barcode tree */
C3SC_barcode_node * srch_C3SCbc_tree(char * str, int i, char * nxt_key, C3SC_barcode_node * scbt_node);

/* set_C3SCbc_tree_revcomp: link C3SC barcode nodes that are the reverse complement of each other */
void set_C3SCbc_tree_revcomp(C3SC_barcode_node * scbt_root);

#endif /* C3SCbc_tree_h */

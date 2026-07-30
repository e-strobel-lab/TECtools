//
//  parse_ct_file.h
//  
//
//  Created by Eric Strobel on 11/20/25.
//

#ifndef parse_ct_file_h
#define parse_ct_file_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../../utils/io_management.h"
#include "../../utils/gen_utils.h"

#include "../../seq_utils/ispair.h"

/* con_table: storage for connectivity table values */
typedef struct con_table {
    char * ttl;   //title
    int len;      //sequence length
    int * n;      //index n
    int * nm1;    //index n-1
    int * np1;    //index n+1
    int * pr2;    //paired nucleotide index n
    int * nat;    //native numbering
    char * bs;    //base identity
    char * db;    //dot-bracket notation
    char * db_an; //dot-bracket notation with base pair type information
    double dG;    //deltaG
    struct con_table * prv; //pointer to previous con_table in doubly linked list
    struct con_table * nxt; //pointer to next con_table in doubly linked list
} con_table;

/* parse_ct_file: manages CT file parsing */
int parse_ct_file(con_table * ct, char * ct_path, char * nm);

/* parse_ct_line1: parse line 1 of CT file */
int parse_ct_line1(con_table * ct, FILE * p_ct, char * nm);

/* parse_ct_body: parse CT file body */
void parse_ct_body(con_table * ct, FILE * p_ct);

/* parse_digit_substring: read a string until next space char, check that all chars are digits, and terminate string */
void parse_digit_substring(char * line, int * i, char * field_nm);

/* parse_char_substring: read a string until next space char, check that all chars are alphabetic, and terminate string */
void parse_char_substring(char * line, int * i, char * field_nm);

/* init_con_table_mem: initialize connectivity table memory */
void init_con_table_mem(con_table * ct);

/* free_con_table_mem: free connectivity table memory */
void free_con_table_mem(con_table * ct, int ct_cnt);

/* ct2dot: convert connectivity table to dot-bracket notation */
void ct2dot(con_table * ct);

#endif /* parse_ct_file_h */

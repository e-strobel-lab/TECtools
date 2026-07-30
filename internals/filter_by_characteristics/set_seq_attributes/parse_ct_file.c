//
//  parse_ct_file.c
//  
//
//  Created by Eric Strobel on 11/20/25.
//

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

#include "parse_ct_file.h"

/* parse_ct_file: manages CT file parsing */
int parse_ct_file(con_table * ct, char * ct_path, char * nm)
{
    FILE * p_ct = NULL;           //file pointer for ct file
    
    //PSEUDOKNOT TEST
    /*
    char test_ttl[5]  = "ex_ct";
    int test_len     = 20;
    int test_n[20]   = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    int test_nm1[20] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    int test_np1[20] = {  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
    int test_pr2[20] = { 14, 13, 12, 11,  0, 20, 19, 18, 17,  0,  4,  3,  2,  1,  0,  0,  9,  8,  7,  6};
    int test_nat[20] = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    
    
    con_table ex_ct = {0};
    ex_ct.ttl = &test_ttl[0];
    ex_ct.len = test_len;
    ex_ct.n   = &test_n[0];
    ex_ct.nm1 = &test_nm1[0];
    ex_ct.np1 = &test_np1[0];
    ex_ct.pr2 = &test_pr2[0];
    ex_ct.dG  = 0.0;
    
    ct = &ex_ct;
    ct2dot(ct);
    printf("%s\n", ct->db);
    abort();
    */
    //END PSEUDOKNOT TEST
    
    int i = 0;                  //general purpose index
    int proceed = 1;            //flag that indicates whether to process another connectivity table within the ct file
    
    con_table * crnt_ct = ct;   //pointer to current con_table structure
    con_table * prev_ct = NULL; //pointer to previous con_table structure
    
    get_file(&p_ct, ct_path); //open ct file
    
    for (i = 0, proceed = 1; proceed == 1; i++) { //until all connectivity tables have been processed
        
        //after the first connectivity table is stored, subsequent connectivity tables
        //from the same ct file are stored as a linked list. memory is allocated for the
        //next connectivity table regardless of whether the ct file contains another
        //connectivity table and then freed if it does not
        if (i >= 1) {
            
            //allocate memory for next connectivity table
            if ((crnt_ct->nxt = malloc(1 * sizeof(*(crnt_ct->nxt)))) == NULL) {
                printf("parse_ct_file: failed to allocate connectivity table memory. aborting...\n");
                abort();
            }
                        
            prev_ct = crnt_ct;      //point prev_ct to the last con_table structure that was processed
            crnt_ct = crnt_ct->nxt; //point crnt_ct to the newly allocated con_table structure
            crnt_ct->prv = prev_ct; //establish doubly linked list
        }
        
        proceed = parse_ct_line1(crnt_ct, p_ct, nm); //parse ct file line 1. the return value indicates whether
                                                     //there is another connectivity table to process
        
        if (proceed) {                    //if ct file contains another connectivity table
            init_con_table_mem(crnt_ct);  //initialize connectivity table memory
            parse_ct_body(crnt_ct, p_ct); //parse ct file body
            ct2dot(crnt_ct);              //store connectivity table as dot-bracket notation
        } else if (i >= 1) {              //if ct file does not contain another connectivity table
            free(prev_ct->nxt);           //free the last con_table struct that was allocated
            prev_ct->nxt = NULL;          //set the pointer to the next con_table to NULL
        }
    }
    
    return i-1; //return number of ct tables parsed from ct file
}

/* parse_ct_line1: parse line 1 of CT file */
int parse_ct_line1(con_table * ct, FILE * p_ct, char * nm)
{
    int i = 0;                      //general purpose index
    char line[MAX_LINE+1] = {0};    //line storage
    char * p_val = NULL;            //pointer to current value string
    char * p_nrg = {0};             //pointer to "ENERGY = " string
    char nrg_str[10] = "ENERGY = "; //energy string for use in comparison
    
    if (get_line(line, p_ct) == 0) { //if there is no line to get, return 0 to indicate that the full ct file was processed
        return 0;
    }
    for (i = 0; isspace(line[i]); i++) {;} //iterate past leading whitespace
    
    //get sequence length
    if (isdigit(line[i])) {                                 //if first char is a digit
        p_val = &line[i];                                   //set value pointer
        parse_digit_substring(line, &i, "sequence length"); //parse seuquence length string
        ct->len = atoi(p_val);                              //store sequence length as an integer
        
    } else { //first seq len char is not a digit, throw error and abort
        printf("parse_ct_file: error - expected first value in ct file first line to be sequence length. aborting...\n");
        abort();
    }
    
    while (isspace(line[i])) {i++;} //iterate past whitespace
    p_val = &line[i];               //set value pointer
    
    //parse delta G value
    p_nrg = strstr(p_val, nrg_str);               //check for energy string
    if ((uint64_t)(p_nrg) == (uint64_t)(p_val)) { //if energy string is at the start of p_val
        p_val = &p_val[strlen(nrg_str)];          //point p_val to the char after energy string
        
        for (i += strlen(nrg_str); !isspace(line[i]) && line[i]; i++) {;} //iterate to next space char
        line[i] = '\0';                                                   //terminate delta G value string
        
        if (check_float_str(p_val, RETURN_OUTCOME)) { //check that delta G value is a float
            ct->dG = strtod(p_val, NULL);             //store delta G value
        } else {
            //TODO: consider adding file name info to error mesages here and elsewhere
            printf("parse_ct_file: error - error when parsing energy value in ct file. aborting...\n");
            abort();
        }
        
        for (i++; isspace(line[i]); i++) {;} //iterate past whitespace
        p_val = &line[i];                    //set pointer to sequence title
        
    } else if (p_nrg == NULL) { //if energy string was not found
        ct->dG = NAN;           //set delta G to NAN
        
    } else { //unexpected format, throw error and abort
        printf("parse_ct_file: error - unexpected format for line 1. aborting...\n");
        abort();
    }
    
    //parse sequence title
    if (!strcmp(p_val, nm)) { //if sequence title matches fasta name, store sequence title
        if ((ct->ttl = malloc((strlen(p_val)+1) * sizeof(*ct->ttl))) == NULL) {
            printf("parse_ct_file: error - failed to allocate memory for ct file title. aborting...\n");
            abort();
        }
        strcpy(ct->ttl, p_val);
        
    } else { //otherwise, throw error and abort
        printf("parse_ct_file: error - sequence title in line 1 of ct file does not match name of input fasta file. aborting...\n");
        abort();
    }
    
    return 1; //return 1 to indicate that there is another connectivity table to process
}

/* parse_ct_body: parse CT file body */
void parse_ct_body(con_table * ct, FILE * p_ct)
{
    int i = 0;                      //general purpose index
    int j = 0;                      //general purpose index
    char line[MAX_LINE+1] = {0};    //line storage
    char * p_val = NULL;            //pointer to current value string
    
    for (i = 0; i < ct->len; i++) { //for every line of the CT file body
        get_line(line, p_ct);       //get body line from connectivity table
        
        for (j = 0; isspace(line[j]); j++) {;} //iterate past leading whitespace
        
        //parse index n
        p_val = &line[j];
        parse_digit_substring(line, &j, "index n");
        ct->n[i] = atoi(p_val);
        while (isspace(line[j])) {j++;}
        
        //parse base identity
        p_val = &line[j];
        parse_char_substring(line, &j, "base identity");
        if (strlen(p_val) == 1) {
            ct->bs[i] = p_val[0];
        } else {
            printf("parse_ct_body: error - expected base identity string in ct file body line to be 1 character. aborting...");
            abort();
        }
        while (isspace(line[j])) {j++;}
        
        //parse index n-1
        p_val = &line[j];
        parse_digit_substring(line, &j, "index n-1");
        ct->nm1[i] = atoi(p_val);
        while (isspace(line[j])) {j++;}
        
        //parse index n+1
        p_val = &line[j];
        parse_digit_substring(line, &j, "index n+1");
        ct->np1[i] = atoi(p_val);
        while (isspace(line[j])) {j++;}
        
        //parse index n of paired nucleotide
        p_val = &line[j];
        parse_digit_substring(line, &j, "index n of paired nucleotide");
        ct->pr2[i] = atoi(p_val);
        while (isspace(line[j])) {j++;}
        
        //parse natural numbering
        p_val = &line[j];
        parse_digit_substring(line, &j, "natural numbering");
        ct->nat[i] = atoi(p_val);
        while (isspace(line[j])) {j++;}
    }
    
    return;
}

/* parse_digit_substring: read a string until next space char, check that all chars are digits, and terminate string */
void parse_digit_substring(char * line, int * i, char * field_nm)
{
    while (!isspace(line[*i]) && line[*i]) { //until whitespace or the string end is reached
        if (isdigit(line[*i])) {             //if the current character is a digit
            (*i)++;                          //increment i
        } else {                             //otherwise, throw error and abort
            printf("parse_digit_substring: expected %s string to be composed of digits. aborting...\n", field_nm);
            
            printf("%c %d\n", line[*i], line[*i]);
            abort();
        }
    }
    
    line[(*i)++] = '\0'; //terminate digit string
    
    return;
}

/* parse_char_substring: read a string until next space char, check that all chars are alphabetic, and terminate string */
void parse_char_substring(char * line, int * i, char * field_nm)
{
    while (!isspace(line[*i]) && line[*i]) { //until whitespace or the string end is reached
        if (isalpha(line[*i])) {             //if the current character is alphabetic
            (*i)++;                          //increment i
        } else {                             //otherwise, throw error and abort
            printf("parse_char_substring: expected %s string to be composed of characters. aborting...\n", field_nm);
            abort();
        }
    }
    
    line[(*i)++] = '\0'; //terminate digit string
    
    return;
}

/* init_con_table_mem: initialize connectivity table memory */
void init_con_table_mem(con_table * ct)
{
    //allocate memory for con_table members
    ct->n   = malloc(ct->len * sizeof(*ct->n));
    ct->nm1 = malloc(ct->len * sizeof(*ct->nm1));
    ct->np1 = malloc(ct->len * sizeof(*ct->np1));
    ct->pr2 = malloc(ct->len * sizeof(*ct->pr2));
    ct->nat = malloc(ct->len * sizeof(*ct->nat));
    ct->bs  = malloc(ct->len+1 * sizeof(*ct->bs));
    
    //if memory allocation for any member failed, throw error and abort
    if (ct->n   == NULL ||
        ct->nm1 == NULL ||
        ct->np1 == NULL ||
        ct->pr2 == NULL ||
        ct->nat == NULL ||
        ct->bs  == NULL ) {
        printf("init_con_table_mem: error - failed to allocate memory for connectivity table contents. aborting...\n");
        abort();
    }
    
    return;
}

/* free_con_table_mem: free connectivity table memory */
void free_con_table_mem(con_table * ct, int ct_cnt)
{
    int i = 0;                  //general purpose index
    con_table * crnt_ct = ct;   //current con_table being freed
    con_table * prev_ct = NULL; //previous con_table
    
    //iterate through the con_table linked list and free all structure members
    
    for (i = 0; i < ct_cnt; i++) {
        free(crnt_ct->ttl);
        free(crnt_ct->n);
        free(crnt_ct->nm1);
        free(crnt_ct->np1);
        free(crnt_ct->pr2);
        free(crnt_ct->nat);
        free(crnt_ct->bs);
        free(crnt_ct->db);
        free(crnt_ct->db_an);
        
        if (i+1 < ct_cnt) {         //if there is another con_table in the linked list
            crnt_ct = crnt_ct->nxt; //point crnt_ct to the next con_table
        }
    }
    
    //free con_table structs that were allocated during linked list generation
    if (crnt_ct->nxt != NULL) { //check that the last con_table in the linked list was reached in the loop above
        printf("free_con_table_mem: error - did not reach the end of the linked list when freeing connectivity table memory. aborting\n");
        abort();
        
    } else if (ct_cnt > 1) {            //if there was more than one connectivity table in the ct file
        while (crnt_ct->prv != NULL) {  //until the root con_table struct is reached
            crnt_ct = crnt_ct->prv;     //set crnt_ct to point to the previous con_table in the linked list
            free(crnt_ct->nxt);         //free the next con_table in the linked list
            crnt_ct->nxt = NULL;        //set the pointer to the next con_table to NULL
        }
    }
    
    return;
}

/* ct2dot: convert connectivity table to dot-bracket notation */
void ct2dot(con_table * ct)
{
    //printf("in ct2dot\n");
    
    int i = 0; //general purpose index
    
    //last_pr2 points to the last value of pr2. this is used to test whether
    //the pr2 value of the current array index is less than or greater than
    //the last pr2 value. if pr2 is less than last_pr2, this closes another
    //base pair in the stem. if pr2 is greater than last_pr2, this indicates
    //the presence of a pseudoknot
    int last_pr2_init = 0;           //initializer for last_pr2 pointer
    int * last_pr2 = &last_pr2_init; //pointer to the most recent pr2 value
    
    //closing_hlx_pr points to pr2 of the first nucleotide in a helix. this
    //is used to determine when all nucleotides that participate in a helix
    //have been read so that a new helix can be started
    int closing_hlx_pr_init = INT_MAX;          //initializer for closing_hlx_pr pointer
    int *closing_hlx_pr = &closing_hlx_pr_init; //pointer to pr2 of the first nt in a helix
    
    //variables associated with pseudoknot annotation. the upstream nucleotides
    //that participate in a pseudoknot are indiciated by a lower case letter;
    //the downstream nucleotides that participate in a pseudoknot are indicated
    //by an upper case letter
    int  pk_ix = 0;                                //index of current letter used
    char pk_bp[27] = "abcdefghijklmnopqrstuvwxyz"; //array of pseudoknot indicators
    char * crnt_pk_bp = &pk_bp[pk_ix];             //current pseudoknot indicator
    int inc_pk_ix = 0;                             //flag to increment pk_ix
    
    if ((ct->db = calloc(ct->len+1, sizeof(ct->db))) == NULL) {
        printf("ct2dot: error - failed to allocate dot-bracket notation memory\n");
        abort();
    }
    
    if ((ct->db_an = calloc(ct->len+1, sizeof(ct->db_an))) == NULL) {
        printf("ct2dot: error - failed to allocate dot-bracket notation memory\n");
        abort();
    }

    int pr_typ = 0; //pair type
    
    //convert connectivity table into dot-bracket notation
    for (i = 0; i < ct->len; i++) { //for every nucleotide
                
        if (ct->pr2[i] == 0) {   //if nucleotide is unpaired,
            ct->db[i] = '.';     //print dot to standard dot-bracket
            ct->db_an[i] = '.';  //print dot to annotated base pair dot-bracket
            
        } else if (!ct->db[i]) { //if nucleotide is paired and db has not been set
                        
            if (!(*last_pr2) || ct->pr2[i] < *last_pr2) { //if initiating or continuing helix
                
                //standard dot-bracket
                ct->db[i] = '(';                       //set upstream bracket
                ct->db[ct->pr2[i]-1] = ')';            //set downstream bracket
                
                //annotated base pair dot-bracket
                pr_typ = ispair(ct->bs[i], ct->bs[ct->pr2[i]-1]);
                switch (pr_typ) {
                    case GC_PAIR:
                        ct->db_an[i] = '[';            //set upstream bracket
                        ct->db_an[ct->pr2[i]-1] = ']'; //set downstream bracket
                        break;
                        
                    case AU_PAIR:
                        ct->db_an[i] = '(';            //set upstream bracket
                        ct->db_an[ct->pr2[i]-1] = ')'; //set downstream bracket
                        break;
                        
                    case GU_PAIR:
                        ct->db_an[i] = '{';            //set upstream bracket
                        ct->db_an[ct->pr2[i]-1] = '}'; //set downstream bracket
                        break;
                        
                    default:
                        printf("ct2dot: error - expected connected bases within connectivity table to be able to form a pair. aborting...\n");
                        abort();
                        break;
                }
                
                
                last_pr2 = &ct->pr2[i];                       //set last_pr2 to pairing partner value
                
                if (*closing_hlx_pr == closing_hlx_pr_init) { //if closing_hlx_pr was not set
                    closing_hlx_pr = &ct->pr2[i];             //point closing_hlx_pr to pr2
                }
                
                if (pk_ix) {                                  //if a pseudoknot was being annotated
                    inc_pk_ix = 1;                            //set flag to increment PK indicator index
                }
                
            } else if (ct->pr2[i] > *last_pr2) {              //pr2 > lastpr2 indicates a pseudoknot
                                
                if (inc_pk_ix) {                              //if a previous pseudoknot was fully set
                    if (pk_ix < strlen(pk_bp)) {              //if unused PK indicators remain
                        crnt_pk_bp = &pk_bp[++pk_ix];         //point crnt_pk_bp to the next indicator
                        inc_pk_ix = 0;                        //reset the inc_pk_ix flag to false
                    } else {
                        printf("ct2dot: error - too many pseudoknots for notation scheme. aborting...\n");
                        abort();
                    }
                }
                  
                //standard dot-bracket
                ct->db[i] = *crnt_pk_bp;                      //set upstream nt to lower case indicator
                ct->db[ct->pr2[i]-1] = toupper(*crnt_pk_bp);  //set downstream nt to upper case indicator
                
                //annotated base pair dot-bracket
                ct->db_an[i] = *crnt_pk_bp;                      //set upstream nt to lower case indicator
                ct->db_an[ct->pr2[i]-1] = toupper(*crnt_pk_bp);  //set downstream nt to upper case indicator
            }
            
        } else if (ct->n[i] == *closing_hlx_pr) {  //if the closing nt of the current helix was reached
            last_pr2 = &last_pr2_init;             //set last_pr2 to initializer
            closing_hlx_pr = &closing_hlx_pr_init; //set closing_hlx_pr to initializer
            
        } else {
            ; //dot-bracket character was already set, skip
        }
    }
    
    return;
}

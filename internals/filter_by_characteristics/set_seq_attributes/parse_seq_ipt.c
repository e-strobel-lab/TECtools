//
//  parse_seq_ipt.c
//  
//
//  Created by Eric Strobel on 11/20/25.
//

#include <stdio.h>
#include <stdlib.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"
#include "../../utils/io_management.h"

#include "../../seq_utils/isDNAbase.h"
#include "../../seq_utils/isRNAbase.h"

#include "../filter_by_characteristics_defs.h"
#include "../filter_by_characteristics_structs.h"

#include "parse_seq_ipt.h"

/* parse_seq_ipt: parse input multiple sequence alignment file */
//note: aligment characters are preseved in the stored sequences and removed later as needed when generating descriptor sequences
int parse_seq_ipt(char * ipt_path, int msa_typ, sequence_attributes ** sq_att) {
        
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int cnt = 0; //number of sequences
    
    char line[MAX_LINE+1] = {0}; //input line storage
    char * nm_ptr = NULL;        //name pointer
    char * sq_ptr = NULL;        //sequence pointer
    
    FILE * ifp = NULL; //input file pointer
    
    //constant strings for line identification
    char var_str[10] = "variants:";
    char wt_str[4]   = "WT:";
    char ref_str[5]  = "REF:";
    
    //open file and count lines
    if ((ifp = fopen(ipt_path, "r")) == NULL) {
        printf("parse_seq_ipt: error - failed to open input sequence file. aborting...\n");
        abort();
    } 
    
    while (get_line(line, ifp)) { //for every line
        if ((msa_typ == STO_FILE && (line[0] != '#' && line[0] != '/')) ||
            (msa_typ == VMT_FILE && (memcmp(line, var_str, strlen(var_str)) &&
                                     memcmp(line, wt_str,  strlen(wt_str))  &&
                                     memcmp(line, ref_str, strlen(ref_str))))) { //if current line is not a header line
            cnt++;                                                               //increment cnt
        }
    }
    
    if (fclose(ifp) == EOF) { //close file
        printf("parse_seq_ipt: error -failed to close input sequence file. aborting...\n");
        abort();
    }
    
    //allocate sequence attribute structure memory
    if ((*sq_att = calloc(cnt, sizeof(**sq_att))) == NULL) {
        printf("parse_seq_ipt: error - failed to allocate memory for sequence attribute structures. aborting...\n");
        abort();
    }
    
    //open input file a second time
    if ((ifp = fopen(ipt_path, "r")) == NULL) {
        printf("parse_seq_ipt: error - failed to open input sequence file. aborting...\n");
        abort();
    }
    
    //get lines, terminate name string, and set pointer to start of sequence string
    i = 0;
    while (get_line(line, ifp) && i < cnt) {
                
        if ((msa_typ == STO_FILE && (line[0] != '#' && line[0] != '/')) ||
            (msa_typ == VMT_FILE && (memcmp(line, var_str, strlen(var_str)) &&
                                     memcmp(line, wt_str,  strlen(wt_str))  &&
                                     memcmp(line, ref_str, strlen(ref_str))))) { //if line is not a header line
            for (j = 0; isspace(line[j]) && line[j] && j < MAX_LINE; j++) {;}    //iterate past leading whitespace
            
            if (line[j]) {         //if loop ended on a char
                nm_ptr = &line[j]; //set name pointer
            } else {               //otherwise, throw error and abort
                printf("parse_seq_ipt: error - unexpected format for sequence line. aborting...\n");
                abort();
            }
            
            //TODO: is this how we want to handle the pipe chars in the msa? Do we want to sub all non-alphanumeric characters with '_'?
            //iterate to end of sequence name, substitute pipe characters with '_'
            while (!isspace(line[j]) && line[j] && j < MAX_LINE) {
                if (line[j] == '|') {
                    line[j] = '_';
                }
                j++;
            }
            
            if (isspace(line[j])) { //if loop ended on space char
                line[j] = '\0';     //terminate name string
            } else {                //otherwise, throw error and abort
                printf("parse_seq_ipt: error - unexpected format for sequence line. aborting...\n");
                abort();
            }
                        
            for (j++; isspace(line[j]) && line[j] && j < MAX_LINE; j++) {;} //iterate past whitespace
            
            if (line[j]) {         //if loop ended on a char
                sq_ptr = &line[j]; //set seq pointer
            } else {               //otherwise, throw error and abort
                printf("parse_seq_ipt: error - line does not contain a sequence");
                abort();
            }
            
            //check that sequence is composed of valid characters
            for (j = 0; sq_ptr[j]; j++) {
                if (!isDNAbase(sq_ptr[j]) && //if char is not a DNA base
                    !isRNAbase(sq_ptr[j]) && //or an RNA base
                    sq_ptr[j] != '.'      && //or a '.' spacer
                    sq_ptr[j] != '-') {      //or a '-' spacer, throw error and abort
                    printf("parse_seq_ipt: error - invalid character '%d' in sequence. aborting...\n", sq_ptr[j]);
                    abort();
                }
            }
                    
            //allocate memory for sequence name
            if (((*sq_att)[i].nm = malloc((strlen(nm_ptr)+1) * sizeof(*((*sq_att)[i].nm)))) == NULL) {
                printf("parse_seq_ipt: error - memory allocation for sequence name failed. aborting...\n");
                abort();
            }
            
            //allocate memory for sequence. len is +2 to account for terminating null and
            //the sequence string starting at 1 so that natural sequence numbering can be
            //used when accessing array members
            if (((*sq_att)[i].sq0 = malloc((strlen(sq_ptr)+2) * sizeof(*((*sq_att)[i].sq0)))) == NULL) {
                printf("parse_seq_ipt: error - memory allocation for sequence name failed. aborting...\n");
                abort();
            }
            
            strcpy((*sq_att)[i].nm, nm_ptr);         //store name in sequence attributes struct
            (*sq_att)[i].sq0[0] = '>';               //set first sq0 array member to '>' as an offset
            (*sq_att)[i].sq1 = &(*sq_att)[i].sq0[1]; //set sq1 pointer to index 1 of sq0
            strcpy((*sq_att)[i].sq1, sq_ptr);        //store sequence starting at index 1 of sq0 (using sq1 ptr)
                        
            //printf("%s\t%s\n", (*sq_att)[i].nm, (*sq_att)[i].sq0);
            //printf("%s\t %s\n\n", (*sq_att)[i].nm, (*sq_att)[i].sq1);
            
            i++; //increment i
        }
    }
    
    //close input file
    if (fclose(ifp) == EOF) { //close file
        printf("parse_seq_ipt: error - failed to close input sequence file. aborting...\n");
        abort();
    }
    
    return cnt; //return sequence count
}

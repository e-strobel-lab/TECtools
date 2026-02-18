//
//  parse_fasta.c
//  
//
//  Created by Eric Strobel on 1/17/23.
//

#include "parse_fasta.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/io_management.h"
#include "./isDNAbase.h"


//TODO: add argument specifying whether the input sequence is expected to be DNA/RNA/contain degenerate base symbols
/* parse_fasta: place name and sequence lines of a fasta file into character arrays */
int parse_fasta(FILE * fp_fasta, char * name, char * seq)
{
    int i = 0;
    int j = 0;
    
    char line[MAX_LINE] = {0}; //array to temporarily store fasta file lines
    
    //get name line from input fasta file
    if (get_line(&line[0], fp_fasta)) {
        if (line[0] == '>') { //fasta line 1 is expected to start with > char
            
            //copy fasta name line to name array
            for (i = 1, j = 0; line[i] && i < MAX_LINE-1; i++, j++) {
                if (!isspace(line[i])) { //spaces are not allowed in the sequence name
                    name[j] = line[i];
                } else { //sequence name contained a space character
                    printf("parse_fasta: error - fasta line contains a space character. aborting...\n");
                    abort();
                }
            }
            name[j] = '\0';
            
            //test that name line doesn't exceed array size - get_line also does this, but doesn't hurt
            //to have it here
            if (i == MAX_LINE-1 && line[i]) {
                printf("parse_fasta: error - fasta file name line exceeded maximum allowed length (%d). aborting...\n", MAX_LINE);
            }
            
        } else { //fasta line 1 did not start with > character, abort
            printf("parse_fasta: error - expected leading '>' in fasta line 1. aborting...\n");
            abort();
        }
    } else { //input fasta file was empty
        return 0;
    }
    
    
    //get sequence line from input fasta file
    if (get_line(&line[0], fp_fasta)) {
        for (i = 0; line[i] && i < MAX_LINE-1; i++) {
            if (isDNAbase(line[i]) || line[i] == 'U' || line[i] == 'u') { //input sequence line should only contain standard DNA char (not allowing whole IUPAC code currently)
                seq[i] = toupper(line[i]);    //set to uppercase when copying line
            } else {
                printf("parse_fasta: error - fasta file sequence line contained the non-standard DNA base %c. the sequence should only contain A/T/G/C/a/t/g/c aborting...\n", line[i]);
                abort();
            }
        }
        seq[i] = '\0';
        
        //test that sequence line doesn't exceed array size - get_line also does this, but doesn't hurt
        //to have it here 
        if (i == MAX_LINE-1 && line[i]) {
            printf("parse_fasta: error - fasta file sequence line exceeded maximum allowed length (%d). aborting...\n", MAX_LINE);
        }
        
    } else { //no second line in input fasta file
        printf("parse_fasta: error - failed to get second fasta file line. aborting...\n");
        abort();
    }
    
    return 1;
}

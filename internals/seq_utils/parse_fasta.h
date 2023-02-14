//
//  parse_fasta.h
//  
//
//  Created by Eric Strobel on 1/17/23.
//

#ifndef parse_fasta_h
#define parse_fasta_h

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "../utils/io_management.h"
#include "./isDNAbase.h"

/* parse_fasta: place name and sequence lines of a fasta file into character arrays
 FILE * fp_fasta: pointer to input fasta file
 char * name: pointer to array for storing sequence name
 char * seq: pointer to array for storing sequence
 */
int parse_fasta(FILE * fp_fasta, char * name, char * seq);

#endif /* parse_fasta_h */

//
//  mk_fasta.c
//  
//
//  Created by Eric Strobel on 1/17/23.
//

#include "mk_fasta.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"
#include "../global/global_structs.h"

#include "./isIUPACbase.h"

char fa_sffx[4] = ".fa";

/* mk_fasta_file: generate fasta file */
int mk_fasta_file(char * nm, char * sq, char * path)
{
    extern char fasta_suffix[4];
    
    FILE * fa_ofp = NULL;         //output fasta file pointer
    char fa_fn[MAX_LINE+1] = {0}; //output fasta filename
    
    int i = 0;   //general purpose index
    int ret = 0; //return value for error checking
    
    for (i = 0; sq[i]; i++) {
        if (!isIUPACbase(sq[i])) {
            printf("mk_fasta_file: error - input sequence contains non-IUPAC DNA base %c. aborting...\n", sq[i]);
            abort();
        }
    }
    
    //generate output fasta file
    if (path != NULL) {
        if (path[strlen(path)-1] != '/') {
            printf("mk_fasta_file: error - last character of path must be a forward slash. aborting...\n");
            abort();
        }
        ret = snprintf(fa_fn, MAX_LINE, "%s%s%s", path, nm, fa_sffx);
    } else {
        ret = snprintf(fa_fn, MAX_LINE, "%s%s", nm, fa_sffx);
    }
    
    if (ret >= MAX_LINE || ret < 0) {
        printf("mk_fasta_file: error - output fasta file path is too long. aborting...\n");
        abort();
    }
    
    //printf("%s\n", fa_fn);
    
    if ((fa_ofp = fopen(fa_fn, "w")) == NULL) {  //open fasta output file
        printf("mk_fasta_file: error - failed to open fasta output file. aborting...");
        abort();
    }
    
    fprintf(fa_ofp, ">%s\n", nm); //print name line to fasta file
    fprintf(fa_ofp, "%s\n", sq);  //print sequence line to fasta file
    
    if ((fclose(fa_ofp)) == EOF) {
        printf("mk_fasta_file: error - error occurred when closing fasta output file. aborting...\n");
        abort();
    }
    
    return 1;
}

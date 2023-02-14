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

/* mk_fasta_file: generate fasta file */
int mk_fasta_file(char * nm, char * sq)
{
    FILE * fa_ofp = NULL;       //output fasta file pointer
    char fa_fn[MAX_LINE] = {0}; //output fasta filename
    
    int i = 0;
    
    for (i = 0; sq[i]; i++) {
        if (!isIUPACbase(sq[i])) {
            printf("mk_fasta_file: error - input sequence contains non-IUPAC DNA base %c. aborting...\n", sq[i]);
            abort();
        }
    }
    
    
    if ((strlen(nm) + 4) >= MAX_LINE) { //check that filename fits in fasta_nm array
        printf("why is your file name so long?\n");
        abort();
    } else {
        //generate output fasta file
        sprintf(fa_fn, "%s.fa", nm);                 //fasta filename
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
    }
    
    return 1;
}

//
//  mk_intermed_trgts.c
//  
//
//  Created by Eric Strobel on 3/17/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "../../global/global_defs.h"
#include "../../global/global_structs.h"

#include "../cotrans_preprocessor_defs.h"
#include "../cotrans_preprocessor_structs.h"

#include "mk_intermed_trgts.h"

/* mk_intermed_trgts: generate individual fasta files for every intermediate
 length of the input sequence, starting at min_len+1 */
int mk_intermed_trgts(char name[MAX_LINE], char seq[MAX_LINE], int min_len)
{
    FILE * out_fp = {NULL};	//FILE pointer for generating output fasta files
    FILE * qc_fp = {NULL};	//FILE pointer for generating QC file
    
    char out_name[MAX_LINE] = {0};	//array for storing output fasta file name
    char qc_name[MAX_LINE] = {0};	//array for storing QC file name
    
    //open file that will contain a list of all targets
    sprintf(qc_name, "./%s_targets/metrics/%s_cotrans_targets_qc.txt", name, name);
    if ((qc_fp = fopen(qc_name,"w")) == NULL) {
        printf("mk_intermed_trgts: error - could not qc output file. Aborting program...\n");
        abort();
    }
    
    static const char ldr[21]  = "atggccttcgggccaa"; //leader sequence, masked by lowercase
    
    int i = 0;
    int j = 0;
    
    int ipt_len = strlen(seq);		//length of input sequence
    char target[MAX_LINE] = {0};	//array to store target
    
    /* starting at min_len+1, generate a fasta file for each
     intermediate transcript length up to the full length
     input sequence
     
     note that i is initialized to 1 because min_len is 1 shorter
     than the shortest target that will be generated */
    for (i = 1; i <= ipt_len-min_len; i++) {
        //copy input sequence as target until
        //current sequence length is reached
        for (j = 0; j < (min_len+i); j++) {
            target[j] = toupper(seq[j]);	//ensure target is uppercase
        }
        target[j] = '\0';
        
        //generate fasta file for current transcript length
        sprintf(out_name, "./%s_targets/intermediate_transcripts/%s_%03dnt_target.fa", name, name, min_len+i);
        if ((out_fp = fopen(out_name,"w")) == NULL) {
            printf("mk_intermed_trgts: error - could not open targets output file. Aborting program...\n");
            abort();
        }
        
        fprintf(out_fp, ">%s_%d_nt\n", name, min_len+i);		//print target name
        fprintf(out_fp, "%s%s\n", ldr, target);					//print leader+target sequence
        
        printf("%03d: %s%s\n",min_len+i, ldr, target);			//print leader+target sequence
        fprintf(qc_fp, "%03d: %s%s\n",min_len+i, ldr, target);	//print leader+target sequence
        
        if ((fclose(out_fp)) == EOF) {
            printf("mk_intermed_trgts: error - error occurred when closing targets output file. Aborting program...\n");
            abort();
        }
    }
    
    if ((fclose(qc_fp)) == EOF) {
        printf("mk_intermed_trgts: error - error occurred when closing qc file. Aborting program...\n");
        abort();
    }
    return 1;
}

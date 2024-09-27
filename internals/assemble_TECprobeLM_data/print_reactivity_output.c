//
//  print_reactivity_output.c
//  
//
//  Created by Eric Strobel on 9/26/24.
//

#include <stdio.h>

#include "../global/global_defs.h"
#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

#include "print_reactivity_output.h"

/* print_reactivity_output: print reactivity values of the enriched transcript lengths*/
void print_reactivity_output(char * out_dir, char * out_nm, mode_parameters * mode_params, int ipt_cnt[TOT_SAMPLES], int nrchd_len[TOT_SAMPLES], double vals[MAX_TRANSCRIPT][TOT_SAMPLES][MAX_IPT], char * seq)
{
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose indep
    
    FILE * out_fp = NULL;        //output file pointer
    char out_fn[MAX_NAME] = {0}; //output file name
    
    //generate reactivity output file
    sprintf(out_fn, "%s/%s_LM_reactivity.txt", out_dir, out_nm); //generate output filename
    if ((out_fp = fopen(out_fn, "w")) == NULL) {                  //open output file
        printf("assemble_TECprobeLM_data: error - failed to open output file. aborting...");
        abort();
    }
    
    //print column headers
    fprintf(out_fp, "nt\tsequence"); //print headers for the nucleotide and sequence columns
        
    //for each input of each sample, print a column header.
    for (i = 0; i < TOT_SAMPLES; i++) {
        for (j = 0; j < ipt_cnt[i]; j++) {
            fprintf(out_fp, "\t%s_S%d_%dnt_%d", out_nm, i+1, nrchd_len[i], j+1);
        }
    }
    fprintf(out_fp, "\n");
        
    //print reactivity values
    for (i = 0; i < nrchd_len[S3] && i < MAX_TRANSCRIPT; i++) {     //for each nucleotide
        fprintf(out_fp, "%d\t%c", i+mode_params->offset, seq[i+1]); //print nucleotide position (index0 is '>')
        for (j = 0; j < TOT_SAMPLES; j++) {                         //for each sample
            for (k = 0; k < ipt_cnt[j]; k++) {                      //for each input
                if (i < nrchd_len[j]) {
                    fprintf(out_fp, "\t%10.6f", vals[i][j][k]);     //print the replicate data value
                } else {
                    fprintf(out_fp, "\t");                          //or print a tab if there is no value
                }
            }
        }
        fprintf(out_fp, "\n");
    }
    
    //close the simple reactivity output file
    if (fclose(out_fp) == EOF) {
        printf("assemble_TECprobeLM_data: error - failed to close output file. Aborting program...\n");
        abort();
    }
}

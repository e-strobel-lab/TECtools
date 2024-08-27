//
//  print_input_filenames.c
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#include "../global/global_defs.h"
#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

#include "print_input_filenames.h"

/* print_input_filenames: print a record of the the input files provided for the current analysis */
void print_input_filenames(char * out_dir, char * out_nm, char * ipt_nm[TOT_SAMPLES][MAX_IPT], int ipt_cnt[TOT_SAMPLES])
{
    int i = 0;                    //general purpose index
    int j = 0;                    //general purpose index
    
    FILE * out_fp = NULL;         //output file pointer
    
    char out_fn[MAX_NAME] = {0};  //outpout file name
    
    //generate output file
    sprintf(out_fn, "%s/%s_inputFiles.txt", out_dir, out_nm); //generate output filename
    if ((out_fp = fopen(out_fn, "w")) == NULL) {              //open output file
        printf("print_input_filenames: error - failed to open output file. aborting...");
        abort();
    }
    
    for (i = 0; i < TOT_SAMPLES; i++) {                            //for each sample
        fprintf(out_fp, "input files for sample %d:\n", i+1);      //print a line indicating the current sample
        for (j = 0; j < ipt_cnt[i]; j++) {                         //for each input file of the current sample
            fprintf(out_fp, "FILE %d:\t%s\n", j+1, ipt_nm[i][j]);  //print the name of the input file
        }
        if (i+1 < TOT_SAMPLES) {                                   //if there is another sample after the current sample
            fprintf(out_fp, "\n");                                 //print a newline
        }
    }
    
    //close the output file
    if (fclose(out_fp) == EOF) {
        printf("print_input_filenames: error - failed to close output file. Aborting program...\n");
        abort();
    }
}

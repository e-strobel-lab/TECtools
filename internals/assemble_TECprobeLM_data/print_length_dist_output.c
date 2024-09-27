//
//  print_length_dist_output.c
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#include "../global/global_defs.h"
#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

#include "print_length_dist_output.h"

/* print_length_dist_output: print a length distribution output file */
void print_length_dist_output(char * out_dir, char * out_nm, mode_parameters * mode_params, int ipt_cnt[TOT_SAMPLES], int max_index, double vals[MAX_TRANSCRIPT][TOT_SAMPLES][MAX_IPT])
{
    double average[MAX_TRANSCRIPT][TOT_SAMPLES] = {0}; //array to store average values for each sample
    
    int i = 0;                    //general purpose index
    int j = 0;                    //general purpose index
    int k = 0;                    //general purpose index
    
    double crrnt_sum = 0;         //sum of the values for the current sample
    
    FILE * out_fp = NULL;         //output file pointer
    char out_fn[MAX_NAME] = {0};  //output file name
    
    //generate output files
    sprintf(out_fn, "%s/%s_LM_lenDist.txt", out_dir, out_nm); //generate output filename
    if ((out_fp = fopen(out_fn, "w")) == NULL) {               //open  output file
        printf("print_length_dist_output: error - failed to open  output file. aborting...");
        abort();
    }
    
    //print headers
    fprintf(out_fp, "transcript");                             //line id is transcrit length
    for (i = 0; i < TOT_SAMPLES; i++) {                        //for each sample
        for (j = 0; j < ipt_cnt[i]; j++) {                     //for each sample input
            fprintf(out_fp, "\t%s_S%d_%d", out_nm, i+1, j+1);  //print a header that contains the sample and replicate number
        }
        fprintf(out_fp, "\t%s_S%d_avg", out_nm, i+1);          //print a header for the current sample average column
    }
    fprintf(out_fp, "\n");                                     //terminate the line with a newline
    
    //print the data lines
    for (i = 0; i < max_index && i < MAX_TRANSCRIPT; i++) {    //while i is less than the max observed line index and MAX_TRANSCRIPT
        fprintf(out_fp, "%d", i+mode_params->offset);          //print the current transcript length
        for (j = 0; j < TOT_SAMPLES; j++) {                    //for each sample
            for (k = 0, crrnt_sum = 0; k < ipt_cnt[j]; k++) {  //for each input of the current sample
                fprintf(out_fp, "\t%10.6f", vals[i][j][k]);    //print the value of the input
                crrnt_sum += vals[i][j][k];                    //add the current input value to the sum of input values for the current sample
            }
            
            average[i][j] = crrnt_sum/ipt_cnt[j];              //calculate the average value of the current sample
            fprintf(out_fp, "\t%10.6f", average[i][j]);        //print the average value
        }
        
        fprintf(out_fp, "\n");                                 //terminate the line with a newline
    }
   
    //close the output file
    if (fclose(out_fp) == EOF) {
        printf("print_length_dist_output: error - failed to close output file. Aborting program...\n");
        abort();
    }
}

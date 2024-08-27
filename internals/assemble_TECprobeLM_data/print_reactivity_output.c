//
//  print_reactivity_output.c
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#include "../global/global_defs.h"
#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

#include "print_reactivity_output.h"

/* print_reactivity_output: print a reactivity output file */
void print_reactivity_output(char * out_dir, char * out_nm, mode_parameters * mode_params, int ipt_cnt[TOT_SAMPLES], int nrchd_len[TOT_SAMPLES], double vals[MAX_TRANSCRIPT][TOT_SAMPLES][MAX_IPT])
{
    double average[MAX_TRANSCRIPT][TOT_SAMPLES] = {0}; //array to store average values for each sample
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose indep
    
    double crrnt_sum = 0; //sum of the values for the current sample inputs
    
    int entry = 0;                             //index of current line group (three lines are printed for each value)
    double LB_offset[3] = {-0.499999, 0, 0.5}; //offset values used to approximate vertical lines
    
    FILE * out_fp = NULL;        //output file pointer
    char out_fn[MAX_NAME] = {0}; //output file name
    
    //calculate the average reactivity of each sample group
    for (i = 0; i < nrchd_len[S3] && i < MAX_TRANSCRIPT; i++) {   //until the enriched transcript length of sample 3 is reached
        for (j = 0; j < TOT_SAMPLES; j++) {                       //for each sample
            if (i < nrchd_len[j]) {                               //if the current index is less than the current enriched length
                for (k = 0, crrnt_sum = 0; k < ipt_cnt[j]; k++) { //sum the reactivity value of each input for the current sample
                    crrnt_sum += vals[i][j][k];
                }
                average[i][j] = crrnt_sum/ipt_cnt[j];             //store the average reactivity
            }
        }
    }
    
    
    //generate output files
    sprintf(out_fn, "%s/%s_LMP_reactivity.txt", out_dir, out_nm); //generate output filename
    if ((out_fp = fopen(out_fn, "w")) == NULL) {                  //open  output file
        printf("print_reactivity_output: error - failed to open output file. aborting...");
        abort();
    }
    
    //print column headers
    fprintf(out_fp, "nt"); //print the header for the nucleotide column
        
    //for each input of each sample, print a column header. then print one
    //additional column header for the average reactivity of the sample
    for (i = 0; i < TOT_SAMPLES; i++) {
        for (j = 0; j < ipt_cnt[i]; j++) {
            fprintf(out_fp, "\t%s_S%d_%dnt_%d", out_nm, i+1, nrchd_len[i], j+1);
        }
        fprintf(out_fp, "\t%s_S%d_%dnt_avg", out_nm, i+1, nrchd_len[i]);
    }
    fprintf(out_fp, "\n");
    
    //print a zero line, which serves as a staring point for the linebars
    fprintf(out_fp, "0.5");                //the linebar begins at 0.5
    for (i = 0; i < TOT_SAMPLES; i++) {    //for each sample
        for (j = 0; j < ipt_cnt[i]; j++) { //for each input
            fprintf(out_fp, "\t0");        //print a zero
        }
        fprintf(out_fp, "\t0");            //then print a zero for the average column
    }
    fprintf(out_fp, "\n");
    
    
    //print data to output file
    for (i = 0; i < nrchd_len[S3] && i < MAX_TRANSCRIPT; i++) { //until the enriched transcript length for sample 3 is reached
        
        //print three entries for each data point:
        //  - entry 1 contains the average reactivity and is the first corner of the bar
        //  - entry 2 contains replicate data (for points) and the average reactivity so the top of the bar is joined in Datagraph
        //  - entry 3 contains the average reactivity and is the second corner of the bar
        
        for (entry = 0; entry < 3; entry++) {                                 //for each line entry
            fprintf(out_fp, "%.6f", i+mode_params->offset+LB_offset[entry]);  //print the line id with the appropriate linebar offset value
            
            for (j = 0; j < TOT_SAMPLES; j++) {                               //for each sample
                if (i < nrchd_len[j]) {                                       //if i is less than the enriched length of the current sample
                    
                    if (entry == 0 || entry == 2) {                           //if printing the linebar corner lines
                        for (k = 0; k < ipt_cnt[j]; k++) {                    //for each input
                            fprintf(out_fp, "\t");                            //print a tab to skip over replicate data values
                        }
                        fprintf(out_fp, "\t%10.6f", average[i][j]);           //then print the average reactivity
                        
                    } else if (entry == 1) {                                  //if printing linebar centers
                        for (k = 0; k < ipt_cnt[j]; k++) {                    //for each input
                            fprintf(out_fp, "\t%10.6f", vals[i][j][k]);       //print the replicate data value
                        }
                        fprintf(out_fp, "\t%10.6f", average[i][j]);           //then print the average reactivity
                        
                    } else {                                                  //sanity check
                        printf("should be impossible to get here\n");
                        abort();
                    }
                    
                } else {                                                      //if the current index is < the current sample enriched length
                    for (k = 0; k < ipt_cnt[j]; k++) {                        //for each input
                        fprintf(out_fp, "\t");                                //print a tab
                    }
                    fprintf(out_fp, "\t");                                    //then print a tab for the average column
                }
            }
            fprintf(out_fp, "\n");                                            //print a newline before proceeding to the next line
        }
    }
   
    //close the output file
    if (fclose(out_fp) == EOF) {
        printf("print_reactivity_output: error - failed to close output file. Aborting program...\n");
        abort();
    }
    
}

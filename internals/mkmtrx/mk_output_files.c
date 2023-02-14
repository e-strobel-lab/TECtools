//
//  mk_output_files.c
//  
//
//  Created by Eric Strobel on 10/11/22.
//

#include "mk_output_files.h"

#include <stdio.h>
#include <stdlib.h>

#include "../global/global_defs.h"

#include "./cotrans_mtrx.h"
#include "./mkmtrx_defs.h"
#include "./mkmtrx_structs.h"

/* print_csv: print csv file of value matrix */
int print_csv(char * prnt_dir, char * smpl_nm, char * output_nm, char * (*mtrx)[MAX_ROW][MAX_COL], int row_min, int row_max, int col_min, int col_max)
{
    FILE * out_fp = NULL; //output reactivity matrix file pointer
    char out_fn[MAX_LINE] = {0}; //output reactivity matrix file name

    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    //generate output files
    sprintf(out_fn, "%s/%s%s", prnt_dir, smpl_nm, output_nm);      //generate reactivity matrix output filename
    if ((out_fp = fopen(out_fn, "w")) == NULL) {                   //open reactivity matrix output file
        printf("print_csv: error - failed to open reactivity output file. aborting...");
        abort();
    }
    
    //print column headers
    fprintf(out_fp, "Transcript Length");

    for (i = col_min; i <= col_max; i++) { //the remaining column headers are nucleotide numbers up to max transcript length
        fprintf(out_fp, ",%d", i);

    }
    fprintf(out_fp, "\n");
    
    for (i = row_min; i <= row_max; i++) { //for every row of the matrix
        //print transcript length number
        fprintf(out_fp, "%d", i);
        
        for (j = col_min; j <= col_max; j++) { //then, for each column
            
            //print the reactivity value if it is present
            if ((*mtrx)[i][j] != NULL) {
                fprintf(out_fp, ",%s", (*mtrx)[i][j]);
            } else {
                fprintf(out_fp, ","); //or print nothing between comma separators if there is no value
            }
        }
        fprintf(out_fp, "\n"); //terminate the row with a newline
    }
    
    //close the output files
    if (fclose(out_fp) == EOF) {
        fprintf(out_fp, "print_csv: error - failed to close output file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

/* columnize_mtrx: transform each matrix into a single column of reactivity values
by concatenating the values for each transcript length */
int columnize_mtrx(char * prnt_dir, char * smpl_nm, cotrans_matrix * mtrx)
{
    int i = 0;          //general purpose index
    int j = 0;          //general purpose index

    FILE * out_fp = NULL; //output file pointer
    char out_fn[MAX_LINE] = {0}; //output file name
    
    sprintf(out_fn, "%s/%s_columns.txt", prnt_dir, smpl_nm); //generate output filename
    if ((out_fp = fopen(out_fn, "w")) == NULL) {             //open output file
        printf("columnize_mtrx: error - failed to open reactivity output file. aborting...");
        abort();
    }
    fprintf(out_fp, "%s\n", smpl_nm);
    
    //TODO: test that enrichment filter is working correctly using a test matrix
    for (i = mtrx->tl[MIN]; i <= mtrx->tl[MAX]; i++) {               //for each transcript
        if (mtrx->Nrchd[i]) {                                        //if the transcript was enriched
            for (j = mtrx->nt[MIN]; mtrx->vals[i][j] != NULL; j++) { //print the reactivity of each nucleotide
                fprintf(out_fp, "%s\n", mtrx->vals[i][j]);
            }
        }
    }
    
    //close the output files
    if (fclose(out_fp) == EOF) {
        fprintf(out_fp, "columnize_mtrx: error - failed to close output file. Aborting program...\n");
        abort();
    }
    return 1;
}


/* print_alignment stats: print summary of alignment states to files and screen.
 per length alignment rates are sent to <smpl_nm>_alignment_rates.txt
 overall aligment totals are setn to <smpl_nm>_alignment_totals.txt */
int print_alignment_stats(char * prnt_dir, char * smpl_nm, alignment_stats * algn, cotrans_matrix * mtrx)
{
    FILE * rates_fp = NULL;        //rates output file pointer
    char rates_fn[MAX_LINE] = {0}; //rates output file name
    
    FILE * ovrll_fp = NULL;        //overall values output file pointer
    char ovrll_fn[MAX_LINE] = {0}; //overall values output file pointer
    
    int i = 0; //general purpose index
    
    /* generate output files */
    sprintf(rates_fn, "%s/%s_alignment_rates.txt", prnt_dir, smpl_nm); //generate output filename string
    if ((rates_fp = fopen(rates_fn, "w")) == NULL) {                   //open output file
        printf("print_alignment_stats: error - failed to open output file. aborting...");
        abort();
    }
    
    sprintf(ovrll_fn, "%s/%s_alignment_totals.txt", prnt_dir, smpl_nm); //generate output filename string
    if ((ovrll_fp = fopen(ovrll_fn, "w")) == NULL) {                    //open output file
        printf("print_alignment_stats: error - failed to open output file. aborting...");
        abort();
    }
    /* end generate output files */
    
    int total_reads[2] = {0};   //variable for counting total reads
    int total_aligned[2] = {0}; //variable for counting aligned reads
    

    for (i = mtrx->tl[MIN]; i <= mtrx->tl[MAX]; i++) { //for every transcript length
        
        total_reads[UNT] += algn[i].tot_C[UNT]; //add total untreated reads for current transcript length to total
        total_reads[MOD] += algn[i].tot_C[MOD]; //add total modified reads for current transcript length to total
        
        total_aligned[UNT] += algn[i].algnd[UNT]; //add aligned untreated reads for current transcript length to total
        total_aligned[MOD] += algn[i].algnd[MOD]; //add aligned modified reads for current transcript length to total
    }
    
    //print header for alignment rates table
    printf("\nlen\t UNTalgnd\t   UNTtot\t UNTrate\t MODalgnd\t   MODtot\tMODrate\t  frac_algnd\n");
    fprintf(rates_fp, "len\t%s_UNTalgnd\t%s_UNTtot\t%s_UNTrate\t%s_MODalgnd\t%s_MODtot\t%s_MODrate\t%s_frac_algnd\n", smpl_nm, smpl_nm, smpl_nm, smpl_nm, smpl_nm, smpl_nm, smpl_nm);

    
    for (i = mtrx->tl[MIN]; i <= mtrx->tl[MAX]; i++) { //for every transcript length
        
        //print alignment rates line for current transcript length
        printf("%d\t%9d\t%9d\t%8.2f\t%9d\t%9d\t%7.2f\t%12.6f\n",
               i,
               algn[i].algnd[UNT], algn[i].tot_C[UNT], algn[i].calc[UNT],
               algn[i].algnd[MOD], algn[i].tot_C[MOD], algn[i].calc[MOD],
               ((float)(algn[i].algnd[UNT] + algn[i].algnd[MOD]))/((float)(total_aligned[UNT] + total_aligned[MOD])));
        fprintf(rates_fp, "%d\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%.6f\n",
               i,
               algn[i].algnd[UNT], algn[i].tot_C[UNT], algn[i].calc[UNT],
               algn[i].algnd[MOD], algn[i].tot_C[MOD], algn[i].calc[MOD],
               ((float)(algn[i].algnd[UNT] + algn[i].algnd[MOD]))/((float)(total_aligned[UNT] + total_aligned[MOD])));


    }
    printf("UNT\t%d\n", total_aligned[UNT]);
    printf("MOD\t%d\n", total_aligned[MOD]);
    printf("%d\n", (total_aligned[UNT] + total_aligned[MOD]));
    
    //print header for aligment totals table
    printf("\n\tuntreated\tmodified\n");
    fprintf(ovrll_fp, "\n\tuntreated\tmodified\n");
    
    //print total reads values to alignment totals table
    printf("reads\t%d\t%d\n", total_reads[UNT], total_reads[MOD]);
    fprintf(ovrll_fp, "reads\t%d\t%d\n", total_reads[UNT], total_reads[MOD]);
    
    //print aligned reads values to alignment totals table
    printf("aligned\t%d\t%d\n\n", total_aligned[UNT], total_aligned[MOD]);
    fprintf(ovrll_fp, "aligned\t%d\t%d\n\n", total_aligned[UNT], total_aligned[MOD]);
    
    //close output files
    if (fclose(rates_fp) == EOF) {
        printf("print_alignment_stats: error - failed to close output file. Aborting program...\n");
        abort();
    }
    
    if (fclose(ovrll_fp) == EOF) {
        printf("print_alignment_stats: error - failed to close output file. Aborting program...\n");
        abort();
    }
    //end close output files
    
    return 1;
}


/* mk_linebars: generate a line that appears as a reactivity bar plot.
 can b used to show overlapping bar plots on the same graph */
int mk_linebars(char * prnt_dir, char * smpl_nm, cotrans_matrix * mtrx)
{
    int i = 0;     //general purpose index
    int j = 0;     //general purpose index
        
    int entry = 0; //index for entering duplicate values in plot, each data point is entered twice
    double offset[2] = {-0.499999,0.5}; //offset values used to approximate vertical lines
    
    FILE * out_fp = NULL;        //output file pointer
    char out_fn[MAX_LINE] = {0}; //output file name
    
    sprintf(out_fn, "%s/%s_linebars.txt", prnt_dir, smpl_nm); //generate output filename
    if ((out_fp = fopen(out_fn, "w")) == NULL) {              //open output file
        printf("mk_linebars: error - failed to open reactivity output file. aborting...");
        abort();
    }
    
    fprintf(out_fp, "nt");                             //print nucleotide column header
    for (i = mtrx->tl[MIN]; i <= mtrx->tl[MAX]; i++) { //print data column headers
        fprintf(out_fp, "\t%s_%d", smpl_nm, i);
    }
    fprintf(out_fp, "\n");
    
    //set value at position 0.5 to zero so that
    //plots start with a vertical line from 0
    fprintf(out_fp, "0.5");
    for (i = mtrx->tl[MIN]; i <= mtrx->tl[MAX]; i++) {
        fprintf(out_fp, "\t0");
    }
    fprintf(out_fp, "\n");
    
    //print data. each value for nt n is printed twice, first
    //at n-0.499999 and second at n+0.5 so that the bar is
    //centerted at position n
    for (i = mtrx->nt[MIN]; i <= mtrx->nt[MAX]; i++) {                //for every transcript
        for (entry = 0; entry < 2; entry++) {                        //loop runs for two offset values
            fprintf(out_fp, "%.6f", (float)(i) + offset[entry]);     //print nt position + offset
            for (j = mtrx->tl[MIN]; j <= mtrx->tl[MAX]; j++) {       //for every transcript
                if (mtrx->vals[j][i] != NULL) {                      //if there's a value
                    fprintf(out_fp, "\t%s", mtrx->vals[j][i]);       //print the value, tab-separated
                } else {
                    fprintf(out_fp, "\t");                           //if no value, print a tab
                }
            }
            fprintf(out_fp, "\n");                                   //end line
        }
    }
    
    //close the output files
    if (fclose(out_fp) == EOF) {
        fprintf(out_fp, "mk_linebars: error - failed to close output file. Aborting program...\n");
        abort();
    }
    return 1;
}

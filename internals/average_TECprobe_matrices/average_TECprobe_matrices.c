//
//  average_TECprobe_matrices.c
//  
//
//  Created by Eric Strobel on 1/16/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "../mkmtrx/cotrans_mtrx.h"

int main(int argc, char *argv[])
{
    FILE * fp_mtrx = NULL;  //matrix file pointer
    
    char output_name[MAX_NAME] = {0}; //ouput file name
    
    cotrans_matrix mtrx = {{0}};         //root cotrans_matrix struct
    cotrans_matrix * crrnt_mtrx = &mtrx; //set crrnt_mtrx pointer to point to root cotrans_matrix
    int mtrx_cnt = 0;                    //counts number of supplied matrices
    
    int i = 0; //general purpose index
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"matrix", required_argument,  0,  'm'},  //matrix file input
            {"output", required_argument,  0,  'o'},  //output file name
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "m:o:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
            case 'm': //matrix supplied
                
                if (mtrx_cnt > 0) { //check if other matrices were supplied previously
                    //if the current matrix file is not the first matrix file, allocate memory
                    //for a new cotrans_matrix struct that extends the cotrans_matrix linked list
                    if ((crrnt_mtrx->nxt = calloc(1, sizeof(*(crrnt_mtrx->nxt)))) == NULL) {
                        printf("main: error - memory allocation for matrix failed. aborting...\n");
                        abort();
                    }
                    crrnt_mtrx = crrnt_mtrx->nxt; //set crrnt_mtrx to point to newly allocated cotrans_mtrx struct
                }
                get_file(&(fp_mtrx), argv[optind-1]);                //set file pointer to matrix file
                get_sample_name(argv[optind-1], &crrnt_mtrx->fn[0]); //set input file name
                strcpy(crrnt_mtrx->nm,crrnt_mtrx->fn);               //set input file name as default matrix name
                store_mtrx(fp_mtrx, crrnt_mtrx);                     //store matrix in cotrans_mtrx struct
                
                //if current matrix is not first matrix, check that row and column
                //counts for the current matrix match those of the first matrix
                //TODO: also check exact numbering?
                if (mtrx_cnt > 0 && (crrnt_mtrx->row_cnt != mtrx.row_cnt || crrnt_mtrx->col_cnt != mtrx.col_cnt)) {
                    printf("main: error - matrix %s does not have the same dimensions as the first matrix (%s). all matrices must have the same dimensions. aborting...\n", crrnt_mtrx->nm, mtrx.nm);
                }
                
                mtrx_cnt++; //increment matrix counter
                break;
                
            case 'o': //output file name supplied
                if (snprintf(output_name, MAX_NAME, "%s.csv", argv[optind-1]) >= MAX_NAME) {
                    printf("average_TECprobe_matrices: error - output name exceeded buffer. aborting...\n");
                    abort();
                }
                break;
            
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    
    if (optind < argc) {
        printf("\nnon-option ARGV-elements:\n");
        while (optind < argc)
            printf("\n%s \n", argv[optind++]);
        putchar('\n');
        printf("Aborting program.\n\n");
        abort();
    }
    /*********** end of option parsing ***********/
    
    if (mtrx_cnt < 2) { //at least two matrix files must be provided
        printf("main: error - no cotranscriptional RNA structure probing matrices were supplied. aborting...\n");
        abort();
    }
    
    if (!output_name[0]) { //no output name was provided, set to default "out.csv"
        strcpy(output_name, "out.csv");
    }
    
    FILE * out_fp = NULL; //output file pointer
    
    //open output file
    if ((out_fp = fopen(output_name, "w")) == NULL) {
        fprintf(out_fp, "average_TECprobe_matrices: error - failed to open output file. Aborting program...\n");
        abort();
    }
    
    int row = 0; //row index
    int col = 0; //column index
    
    double tot = 0; //variable for storing total reactivity value when averaging
        
    fprintf(out_fp, "Transcript Length"); //print transcript length header
    
    //print nucleotide column headers
    for (col = 1; mtrx.vals[0][col] != NULL; col++) {
        fprintf(out_fp, ",%s", mtrx.vals[0][col]);
    }
    fprintf(out_fp, "\n");
    
    
    //calculate average reactivity and print to file
    for (row = 1; row < mtrx.row_cnt; row++) {    //for every row of the matrix

        fprintf(out_fp, "%s", mtrx.vals[row][0]); //print transcript length
        
        for (col = 1; mtrx.vals[row][col] != NULL; col++) { //for each column
            if (mtrx.vals[row][col][0] == '\0') {           //if the matrix position does not contain a value
                fprintf(out_fp, ",");                       //print a ','
                
            } else {                                        //if the matrix position does contain a value
                crrnt_mtrx = &mtrx;                         //set the crrnt_mtrx pointer to the root matrix
                tot = 0;                                    //zero the reactivity total
                
                for (i = 0; i < mtrx_cnt; i++) {                     //for each input matrix
                    tot += strtof(crrnt_mtrx->vals[row][col], NULL); //add the reactivity value of the current matrix position to the total
                    crrnt_mtrx = crrnt_mtrx->nxt;                    //point crrnt_mtrx to the next matrix
                }
                
                fprintf(out_fp, ",%.6f", tot/((float)mtrx_cnt)); //print the average reactivity value
            }
        }
        
        fprintf(out_fp, "\n"); //terminate row with a newline
    }
    
    //close output file
    if (fclose(out_fp) == EOF) {
        printf("average_TECprobe_matrices: error - failed to close output file. Aborting program...\n");
        abort();
    }
}





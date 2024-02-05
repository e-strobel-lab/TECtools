//
//  cotrans_mtrx.c
//  
//
//  Created by Eric Strobel on 1/16/24.
//

#include <stdio.h>

#include "cotrans_mtrx.h"

/* store_mtrx: store cotranscriptional RNA structure probing matrix
 as a multi-dimensional array in cotrans_matrix struct. */
int store_mtrx(FILE * ifp, cotrans_matrix * mtrx)
{
    int i = 0;        //general purpose index
    int j = 0;        //general purpose index
    int k = 0;        //general purpose index
    int line_cnt = 0; //counts number of lines in matrix file
    int more_vals = 1; //flag that matrix line contains another value
    int prev_nt = 0;  //previous nucleotide value in column header
    int crrnt_nt = 0; //current  nucleotide value in column header
    
    int prev_len = 0;  //number of values in previous row of matrix
    int crrnt_len = 0; //expected number of values in currnt row of matrix
    int val_cnt = 0;   //counts number of values in current matrix row
    
    char line[MAX_LINE] = {0};     //array to store matrix line
    char tmp_val[MAX_FIELD] = {0}; //array to store value from matrix line
    
    for (line_cnt = 0; get_line(line, ifp); line_cnt++) { //for every line of the matrix
        
        //matrix is stored as a csv file. the loops below
        //iterate through comma-separated values and copies
        //them into a temporary array while checking that
        //each character within the matrix is valid. memory
        //is allocated for the value within the cotrans_matrix
        //struct, and the string stored in the temporary array
        //is copied to the allocated memory.
        
        //for every value in the data line. more_vals indicates
        //that there is another value in the line
        for (i = 0, j = 0, val_cnt = 0, more_vals = 1; more_vals; i++) {
                                                          
            /* check each character of the value for validity as it is stored in tmp_array.
            tmp_val is zeroed at the start of each loop so that the number of columns that
            contain a value can be tracked.
             */
            for (k = 0, tmp_val[0] = '\0'; line[j] != ',' && line[j] != '\r' && line[j]; j++, k++) {
                if (isdigit(line[j]) || //char must be a digit
                    line[j] == '.'   || //or '.'
                    line[j] == '-'   || //or '-'
                    line[j] == 'e'   || //or 'e'
                    line[j] == 'E') {   //or 'E'
                    tmp_val[k] = line[j];
                    
                } else if (line_cnt && i) {
                    //coordinate 0,0 in matrix can contain any character. all other cells
                    //must only contain the valid characters listed above. throw error and
                    //abort if invalid character is detected
                    printf("store_mtrx: error - unexpected character %c (%d) in input matrix data line. aborting...\n", line[j], line[j]);
                    abort();
                }
            }
            tmp_val[k] = '\0'; //set terminating null char in tmp_val string
            
            
            /* when reading the column header line (row 0), the tests below are performed to confirm
            that nucleotide numbering starts at 1 and that the column headers are numbered sequentially.
            */
            if (!line_cnt && i >= 1) {    //the tests below are not performed for column 0
                crrnt_nt = atoi(tmp_val); //store integer value of tmp_val string
                
                if (i == 1 && crrnt_nt != 1) { //column numbering must start at 1
                    printf("store_mtrx: error - nucleotide numbering in matrix %s starts at %d. nucleotide numbering must start at 1. aborting...\n", mtrx->fn, crrnt_nt);
                    abort();

                } else if (i >= 2 && crrnt_nt != (prev_nt+1)) { //column numbering must be sequential
                    printf("\n%d\t%d\t%d\n", i, crrnt_nt, prev_nt+1);
                    
                    printf("store_mtrx: error - nucleotide numbers in column headers must be sequential. aborting...\n");
                    abort();
                }
                
                prev_nt = crrnt_nt; //set prev col header value for comparison in next loop iteration
            }
            
            
            /* when reading data lines, the value in column 0 is the expected number of data values in
             the subsequent columns of the row and is stored as crrnt_len. for all data lines after the
             first data line, the current row length compared to the previous row length to confirm that
             row number is sequential. note that the last operation of the current loop iteration is to
             set the prev_len variable as crrnt_len for the next loop iteration, since crrnt_len is
             needed for several other tests below
             
             for all subsequent columns, if data is present, val_cnt is incremented to track the total
             number of values for the line.
             */
            if (line_cnt && !i) {           //reading first value in row
                crrnt_len = atoi(tmp_val);  //set expected length of current row
                
                //test that current expected row length is prev_len+1. This
                //is only done after reaching the second data line so that a
                //comparison can be made
                if (line_cnt >= 2 && crrnt_len != (prev_len+1)) {
                    printf("store_mtrx: error - current row length (%d) is not equal to the previous row length (%d) plus one. aborting...\n", crrnt_len, prev_len);
                    abort();
                }
                
            } else if (tmp_val[0]) { //data is present in current field
                val_cnt++;           //increment number of values in current row
            }

            
            /* use length of tmp_val to allocate memory in the corresponding entry
             of the mtrx->vals array. if successful, copy tmp_val string to the newly
             allocated memory, othwerwise throw error and abort
             */
            if ((mtrx->vals[line_cnt][i] = malloc((strlen(tmp_val)+1) * sizeof(mtrx->vals[line_cnt][i]))) == NULL) {
                printf("store_mtrx: error - memory allocation for matrix field value failed. aborting...\n");
                abort();
            } else {
                strcpy(mtrx->vals[line_cnt][i], tmp_val); //copy tmp_val to mtrx-vals array
            }
            
            if (!line[j] || line[j] == '\r') { //reached terminating null or carriage return in line
                more_vals = 0; //set more_vals to zero to exit the loop
            } else {
                j++; //more values in current line, incrment past delimiter to next char
            }
        }
        
        
        /* if current line is column header line, confirm that minimum column count is met and
         store column count, minimum nt value, and maximum nt value in cotrans_matrix struct
         */
        if (!line_cnt) { //reading column header line
            if (i >= MIN_COL) {    //test that number of columns meets minimum requirement
                mtrx->col_cnt = i; //set column count value in cotrans_matrix struct
                mtrx->nt[MIN] = atoi(mtrx->vals[line_cnt][1]);   //set minimum nucleotide value
                mtrx->nt[MAX] = atoi(mtrx->vals[line_cnt][i-1]); //set maximum nucleotide value
            } else {
                printf("store_mtrx: error - expected cotranscriptional RNA structure probing matrix to contain at least %d columns. aborting...\n", MIN_COL);
                abort();
            }
        
        /* if current line is a data line: if at row 0, set minimum transcript length value. for
         all rows, check that the number of columns matches the expected value (determined using
         the column header line) and that the number of data values in the current row matches
         the expected number that is specified by column 0 of the row
         */
        } else if (line_cnt > 0) { //reading data line
            if (line_cnt == 1) {   //data line is first transcript length
                mtrx->tl[MIN] = atoi(mtrx->vals[line_cnt][0]); //set minimum transcript length value
            }
            
            if (i != mtrx->col_cnt) { //check that number of cols matches expected value
                printf("store_mtrx: error - number of fields in data line (%d) does not match the number of column headers (%d). aborting...\n", i, mtrx->col_cnt);
                abort();
            }
            
            if (val_cnt != crrnt_len) { //check that expected number of values was observed for current row
                printf("store_mtrx: error - number of values in current row (%d) is not equal to the expected number (%d) aborting...\n", val_cnt, crrnt_len);
                abort();
            }
            
        } else { //this error should not be possible
            printf("store_mtrx: error - matrix line count is negative. this should never happen. aborting...\n");
            abort();
        }
        
        prev_len = crrnt_len; //set previous row length to current row length before next iteration of loop
    }
    
    
    /* after processing all rows, set the maximum transcript length value to the value of
     column 0 in the final row. then check that the number of rows meets the minimum required
     number of rows and set the row_count variable in the cotrans_matrix struct
     */
    mtrx->tl[MAX] = atoi(mtrx->vals[line_cnt-1][0]); //set maximum transcript length value

    if (line_cnt >= MIN_ROW) {    //test that number of rows meets minimum requirement
        mtrx->row_cnt = line_cnt; //set row count value in cotrans_matrix struct
    } else {
        printf("store_mtrx: error - expected cotranscriptional RNA structure probing matrix to contain at least %d rows. aborting...\n", MIN_ROW);
        abort();
    }
    
    return 1;
}

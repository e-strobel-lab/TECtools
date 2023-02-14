//
//  mtrx2cols.c
//  
//
//  Created by Eric Strobel on 8/26/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>

#include "../global/global_defs.h"
#include "../utils/io_management.h"
#include "../mkmtrx/cotrans_mtrx.h"

#define MAX_FIELD 64     //max array size for storing value from csv file
#define MIN_ROW 2
#define MIN_COL 2
#define MAX_WINDOW 8     //max number of windows that can be specified in filter file

/* structure declarations */
typedef struct filter {             //constraints for filtering cotranscriptional RNA structure probing matrices
    int ncnt;                       //number of nucleotides included in filter
    int wcnt;                       //number of transcript length windows specified by filter
    int nts[MAX_COL];               //array of nucleotides included in filter
    int wndw[MAX_WINDOW][2];        //array of transcript length windows, each window is specified by a min and max value
} filter;

int store_mtrx(FILE * ifp, cotrans_matrix * mtrx);     //read input matrix file and store data in cotrans_mtrx array
int apply_aliases(FILE * ifp, cotrans_matrix * mtrx);  //replace default matrix names with user-provided aliases
int parse_filter(FILE * ifp, filter * fltr);           //parse filter file for included nucleotides and windows
int filter_mtrx(cotrans_matrix * mtrx, filter * fltr, char * output_name); //filter matrices and print output to file

int main(int argc, char *argv[])
{
    FILE * fp_mtrx = NULL;  //matrix file pointer
    FILE * fp_alias = NULL; //alias file pointer
    FILE * fp_fltr = NULL;  //filter file pointer
    
    char output_name[MAX_NAME] = {0}; //ouput file name
    
    cotrans_matrix mtrx = {{0}};         //root cotrans_matrix struct
    cotrans_matrix * crrnt_mtrx = &mtrx; //set crrnt_mtrx pointer to point to root cotrans_matrix
    int mtrx_cnt = 0;                    //counts number of supplied matrices
    int alias_supplied = 0;              //flag that alias file was supplied
    
    filter fltr = {0};      //struct to store filter constraints
    int fltr_supplied = 0;  //flag that filter file was supplied
    
    int i = 0; //general purpose index
    
    
    /****** parse options using getopt_long ******/
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"matrix",        required_argument,  0,  'm'},  //matrix file input
            {"filter",        required_argument,  0,  'f'},  //filter file input
            {"alias",         required_argument,  0,  'a'},  //alias  file input
            {"output",        required_argument,  0,  'o'},  //output file name
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "m:f:a:o:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
                
            case 'm': //matrix supplied
                if (alias_supplied) { //check if alias file was previously supplied
                    //alias file must be supplied only after all matrix
                    //files have been supplied. throw error and abort.
                    printf("main: error - the alias file option must be supplied after all matrix input files. aborting...\n");
                    abort();
                }
                
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
                if (mtrx_cnt > 0 && (crrnt_mtrx->row_cnt != mtrx.row_cnt || crrnt_mtrx->col_cnt != mtrx.col_cnt)) {
                    printf("main: error - matrix %s does not have the same dimensions as the first matrix (%s). all matrices must have the same dimensions. aborting...\n", crrnt_mtrx->nm, mtrx.nm);
                }
                
                mtrx_cnt++; //increment matrix counter
                break;
            
            case 'f': //filter supplied
                get_file(&(fp_fltr), argv[optind-1]);  //set file pointer to filter file
                parse_filter(fp_fltr, &fltr);          //parse filter file and store constraints in fltr struct
                fltr_supplied++;                       //increment number of supplied filter files
                break;
                
            case 'a': //aliases supplied
                get_file(&(fp_alias), argv[optind-1]); //set file pointer to alias file
                apply_aliases(fp_alias, &mtrx);        //replace default matrix names with aliases
                alias_supplied++;                      //increment number of supplied alias files
                break;
                
            case 'o': //output file name supplied
                sprintf(output_name, "%s.txt", argv[optind-1]);
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
    
    if (!mtrx_cnt) { //at least 1 matrix file must be provided
        printf("main: error - no cotranscriptional RNA structure probing matrices were supplied. aborting...\n");
        abort();
    }
    
    if (fltr_supplied != 1) { //1 and only 1 filter file must be provided
        printf("main: error - %d filter files were supplied. expected 1. aborting...\n", fltr_supplied);
        abort();
    }
    
    if (alias_supplied > 1) { //no more than one alias file can be provided
        printf("main: error - %d alias files were supplied. expected a maximum of 1. aborting...\n", alias_supplied);
        abort();
    }
    
    if (!output_name[0]) { //no output name was provided, set to default "out.txt"
        strcpy(output_name, "out.txt");
    }
    
    /* print input details */
    printf("\nprocessing %d matrices:\n", mtrx_cnt);
    crrnt_mtrx = &mtrx; //set crrnt_mtrx to point to root cotrans_mtrx struct
    while (crrnt_mtrx != NULL) {
        printf("  %s\n", crrnt_mtrx->nm);
        crrnt_mtrx = crrnt_mtrx->nxt;
    }
    
    printf("with the following filter:\n");
    printf("  nucleotides:\n   ");
    for (i = 0; i < fltr.ncnt; i++) {
        printf(" %d", fltr.nts[i]);
    }
    printf("\n  in transcript length windows:\n");
    for (i = 0; i < fltr.wcnt; i++) {
        printf("    min=%d,max=%d\n", fltr.wndw[i][MIN], fltr.wndw[i][MAX]);
    }
    /* end print input details */
    
    filter_mtrx(&mtrx, &fltr, output_name); //filter matrices and generate output file
}

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
            for (k = 0, tmp_val[0] = '\0'; line[j] != ',' && line[j]; j++, k++) {
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
            
            if (!line[j]) {    //reached terminating null in line
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
  
/* apply_aliases: set nm variable in cotrans_matrix struct
 to an alias specified by a user-provided alias file*/
int apply_aliases(FILE * ifp, cotrans_matrix * mtrx)
{
    int i = 0; //general purpose index

    cotrans_matrix * crrnt_mtrx = NULL; //pointer to current cotrans_matrix struct
    
    char line[MAX_LINE] = {0};  //array to store aliase file line
    int line_cnt = 0;           //counts number of lines in alias file
    char * p_ipt_nm = NULL;     //pointer to input name string in alias file line
    char * p_alias = NULL;      //pointer to alias string in alias file line
    
    int crrnt_mtch = 0;
    int mtch_cnt = 0;  //tracks number of default matrix name matches
    
    for (line_cnt = 0; get_line(line, ifp); line_cnt++) {
        
        //iterate to first tab in line. check that the preceding
        //string does not contain any space characters
        for (i = 0; line[i] != '\t' && line[i]; i++) {
            if (isspace(line[i])) { //input name cannot contain space characters, throw error and abort
                printf("apply_aliases: error - input name string cannot contain space characters. aborting...");
                abort();
            }
        }
        
        if (line[i] == '\t') {
            //found tab. split line string at a tab and set pointers
            //to the resulting input name and alias strings
            line[i] = '\0';      //set tab to null character
            p_ipt_nm = &line[0]; //set input name pointer to start of line array
            if (line[i+1] && !isspace(line[i+1])) { //if tab is followed by a non-space, non-null character
                p_alias = &line[i+1];               //set alias pointer to index that follows the tab (now null char)
                for (i = 0; p_alias[i]; i++) {      //check that alias string does not contain space characters
                    if (isspace(p_alias[i])) {      //if alias name contains a space character, throw error and abort
                        printf("apply_aliases: error - alias string cannot contain space characters. aborting...\n");
                        abort();
                    }
                }
                if (strlen(p_alias) > MAX_NAME-1) { //check that alias will fit within name array bounds
                    printf("apply_aliases: error - alias string is %lu characters long. Maximum length is %d characters. aborting...\n", strlen(p_alias), MAX_NAME-1);
                    abort();
                }
            } else {
                //alias string does not begin with a non-space, non-null character
                printf("apply_aliases: error - alias string must begin with a non-space, non-null character. aborting...\n");
                abort();
            }
        } else { //no tab was detected, throw error and abort
            printf("apply_aliases: error - expected alias file line to contain a tab separating input name and alias. aborting...");
            abort();
        }
        
        
        //search input matrix linked list for default names that match
        //the alias file input name. if a match is found, replace default
        //name with the alias.
        crrnt_mtrx = mtrx; //set crrnt_mtrx pointer to point to root matrix
        crrnt_mtch = 0;
        while (crrnt_mtrx != NULL) {    //run loop until all matrices are tested
            if (!crrnt_mtrx->alias) {  //alias has not yet been applied to current matrix
                if (!strcmp(crrnt_mtrx->nm, p_ipt_nm)) { //check if alias file input name is a match to default name
                    strcpy(crrnt_mtrx->nm, p_alias);     //if match, replace default name with alias,
                    crrnt_mtrx->alias = 1;               //set flag that alias has been applied
                    crrnt_mtch++;                        //increment number of matches observed for current name
                    mtch_cnt++;                          //increment match counter
                }
            } else if (!strcmp(crrnt_mtrx->nm, p_alias)) {
                printf("apply_aliases: error - alias file contains duplicate entry for alias %s. aborting...\n", p_alias);
                abort();
            }
            crrnt_mtrx = crrnt_mtrx->nxt; //otherwise, proceed to next matrix in the linked list
        }
        
        if (!crrnt_mtch) { //no match was found
            printf("apply_aliases: error - no match for input name %s. aborting...\n", p_ipt_nm);
            abort();
        }
        
        if (crrnt_mtch > 1) { //more than one match was found
            printf("apply_aliases: error - more than one match was found for input name %s. aborting..\n", p_ipt_nm);
            abort();
        }
    }
    
    return 1;
}
 
/* parse_filter: parse filter file for list of nts and transcript length windows
 to use when filtering cotranscriptional RNA structure probing matrices */
int parse_filter(FILE * ifp, filter * fltr)
{
    int i = 0;        //general purpose index
    int j = 0;        //general purpose index
    int k = 0;        //general purpose index
    int more_vals = 1; //flag that nucleotide list contains another value

    char line[MAX_LINE] = {0};     //array to store input line
    char tmp_val[MAX_FIELD] = {0}; //array to temporarily store input line substrings
    
    int wndw_cnt = 0; //counts number of transcript length windows specified by filter file
    
    int nt_seen[MAX_COL] = {0};
    
    /* parse first line of file for nucleotides to include when filtering
     cotranscriptional RNA structure probing matrices */
    
    char nts_id[5] = {"nts="}; //string that indicates line is nts line
    get_line(line, ifp);       //get first line of filter file, which specifies what nucleotides to include
    if (!memcmp(line, nts_id, strlen(nts_id))) { //check for "nts=" line id
        
        j = strlen(nts_id); //initialize j to length of nts id to skip to first value
        
        //copy nts line values to tmp_val array while checking that all non-comma
        //characters are digits. then store tmp_val as an integer in filter struct
        //nts array
        for (i = 0, more_vals = 1; more_vals; i++) { //for each value in the nts line
            
            //iteate through value until ',' char is reached
            for (k = 0, tmp_val[k] = '\0'; line[j] != ',' && line[j]; j++, k++) {
                if (isdigit(line[j])) {   //if char is digit
                    tmp_val[k] = line[j]; //copy char to tmp_val array
                } else { //only digits are allowed in nts values, throw error and abort
                    printf("parse_filter: error - unexpected character in filter file nucleotides line. aborting...\n");
                    abort();
                }
            }
            tmp_val[k] = '\0'; //set terminating null char in tmp_val array
            fltr->nts[i] = atoi(tmp_val); //convert tmp_val string to int, store in filter struct
            if (!nt_seen[fltr->nts[i]]) {  //check if nucleotide was provided more than once in the filter
                nt_seen[fltr->nts[i]] = 1;
            } else {
                printf("parse_filter: error nucleotide %d was included in the filter more than once. aborting...\n", fltr->nts[i]);
                abort();
            }
            
            if (!fltr->nts[i]) {
                printf("parse_filter: error - nucleotide filter must have a value greater than 0. aborting...\n");
                abort();
            }
            
            if (!line[j]) {   //reached end of line string
                more_vals = 0; //set more_vals to zero to exit loop
            } else {
                j++; //line contains another value, increment j index past the comma separator
            }
        }
        
        fltr->ncnt = i; //set nucleotide count in filter struct
        
    } else { //filter file must begin with nucleotides line, throw error and abort
        printf("parse_filter: error - incorrect format for filter file nucleotides line. aborting...\n");
        abort();
    }
    
    /* parse remaining lines in filter file for transcript length windows
     to apply when filtering cotranscriptional RNA structure probing matrices */
    
    char * crrnt_fld = NULL;               //pointer to start of min/max fields in window line
    char wndw_ids[2][5] = {"min=","max="}; //window ids for identifying min and max values
    int id = 0; //index to track number of loop iterations when processing min/max fields
    i = j = k = 0;
    
    
    
    for (wndw_cnt = 0; get_line(line, ifp); wndw_cnt++) { //for each remaining line in filter file
        
        //process line to identify transcript length window min and max values
        for (i = 0, id = 0, crrnt_fld = &line[0]; id < 2 && crrnt_fld[i]; id++) {
            //check that window values are preceded by "min=" and "max=" strings in first,
            //and second iterations of the loop, respectively
            if (!memcmp(crrnt_fld, wndw_ids[id], strlen(wndw_ids[id]))) { //verify that window id is present
                crrnt_fld = &crrnt_fld[strlen(wndw_ids[id])];             //set pointer to start of value
                for (i = 0; crrnt_fld[i] != ',' && crrnt_fld[i]; i++) {
                    if (isdigit(crrnt_fld[i])) {   //if char is digit
                        tmp_val[i] = crrnt_fld[i]; //copy char to tmp_val array
                    } else { //only digits allowed in window values, throw error and abort
                        printf("parse_filter: error - unexpected character in filter file window line. aborting...\n");
                        abort();
                    }
                }
                tmp_val[i] = '\0'; //set teminating null character
                                
                fltr->wndw[wndw_cnt][id] = atoi(tmp_val); //convert tmp_val string to int, store in filter struct
                
                if ((crrnt_fld[i] == ',' && crrnt_fld[i+1])) { //if reached comma separator and next char is not null
                    crrnt_fld = &crrnt_fld[i+1];               //set crrnt_fld to point to next comma separated value
                } else if (!crrnt_fld[i]) {    //if reached end of line
                    crrnt_fld = &crrnt_fld[i]; //set crrnt_fld to point to terminating null
                } else { //line is incorrectly formatted, throw error and abort
                    printf("parse_filter: error - unexpected format for filter window. aborting...");
                    abort();
                }
            }
        }
        
        //check that transcript length window minimum is < maximum
        if (fltr->wndw[wndw_cnt][MIN] > fltr->wndw[wndw_cnt][MAX]) {
            printf("parse_filter: error - transcript length window minimum is greater than transcript window length maximum. aborting...\n");
            abort();
        }
        
    }
    
    fltr->wcnt = wndw_cnt; //store window count in filter struct
    
    return 1;
}

/* filter_mtrx: filter cotranscriptional RNA structure probing matrix using
 user supplied constraint. print values that pass filter to an output file */

int filter_mtrx(cotrans_matrix * mtrx, filter * fltr, char * output_name) {
    
    int m = 0;   //matrix index
    int w = 0;   //window index
    int row = 0; //row index
    int k = 0;   //nucleotide index
    
    FILE * out_fp = NULL; //output file pointer
    
    //open output file
    if ((out_fp = fopen(output_name, "w")) == NULL) {
        fprintf(out_fp, "filter_mtrx: error - failed to open output file. Aborting program...\n");
        abort();
    }
    
    cotrans_matrix * crrnt_mtrx = NULL; //pointer to current cotrans_matrix struct
    
    /* print column headers */
    fprintf(out_fp, "transcript_length");             //print column 0 header
    crrnt_mtrx = mtrx;                                //set pointer to root cotrans_matrix
    while (crrnt_mtrx != NULL) {                      //for every matrix...
        for (w = 0; w < fltr->wcnt; w++) {            //for every window...
            for (k = 0; k < fltr->ncnt; k++) {        //for every nucleotide...
                fprintf(out_fp, "\t%s_%dto%d_nt%d",   //print tab-separated column header that contains
                       crrnt_mtrx->nm,                //matrix name and...
                       fltr->wndw[w][MIN],            //minimum transcript length and...
                       fltr->wndw[w][MAX],            //maximum transcript length and...
                       fltr->nts[k]);                 //nucleotide numberand...
            }
            if (w+1 < fltr->wcnt) {                   //if there are more windows
                fprintf(out_fp, "\t");                //print a tab
            }
        }
        if ((crrnt_mtrx = crrnt_mtrx->nxt) != NULL) { //if there's another matrix
            fprintf(out_fp, "\t");                    //print a tab
        } else {
            fprintf(out_fp, "\n");                    //otherwise print a newline
        }
    }
    
    
    /* print data lines*/
    
    for (crrnt_mtrx = mtrx, row = 1; row < crrnt_mtrx->row_cnt; row++) { //for every row of the matrix
                                                                         //row initializes to 1 to skip column headers
        for (m = 0; crrnt_mtrx != NULL; m++) {                           //for every matrix...
            if (!m) {                                                    //if first matrix,
                fprintf(out_fp, "%s\t", crrnt_mtrx->vals[row][0]);       //print transcript length to column 0
            }
            for (w = 0; w < fltr->wcnt; w++) {                           //for every window...
                for (k = 0; k < fltr->ncnt; k++) {                       //for every nucleotide...
                    
                    //if transcript length is within transcript length window, print value
                    if (atoi(crrnt_mtrx->vals[row][0]) >= fltr->wndw[w][MIN] && //transcript length is >= window minimum
                        atoi(crrnt_mtrx->vals[row][0]) <= fltr->wndw[w][MAX]) { //transcript length is <= window maximum
                        fprintf(out_fp, "%s", crrnt_mtrx->vals[row][fltr->nts[k]]);
                    }
                    
                    if (k+1 < fltr->ncnt || w+1 < fltr->wcnt) { //if there's another nucleotide or window
                        fprintf(out_fp, "\t");                  //print a tab
                    }
                }
                
                if (w+1 < fltr->wcnt) {    //if there's another window
                    fprintf(out_fp, "\t"); //print another tab
                }
            }
            
            if ((crrnt_mtrx = crrnt_mtrx->nxt) != NULL) { //if there's another matrix
                fprintf(out_fp, "\t\t");                  //print two tabs
            } else {
                fprintf(out_fp, "\n");                    //otherwise print a newline
            }
        }
        
        crrnt_mtrx = mtrx; //set crrnt matrix to matrix root before processing next row
    }
    
    //close output file
    if (fclose(out_fp) == EOF) {
        fprintf(out_fp, "filter_mtrx: error - failed to close output file. Aborting program...\n");
        abort();
    }

    return 1;
}

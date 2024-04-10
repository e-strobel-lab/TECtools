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

#define MAX_WINDOW 8     //max number of windows that can be specified in filter file

/* structure declarations */
typedef struct filter {             //constraints for filtering cotranscriptional RNA structure probing matrices
    int ncnt;                       //number of nucleotides included in filter
    int wcnt;                       //number of transcript length windows specified by filter
    int nts[MAX_COL];               //array of nucleotides included in filter
    int wndw[MAX_WINDOW][2];        //array of transcript length windows, each window is specified by a min and max value
} filter;

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
    
    int ret = 0; //variable for storing snprintf return value
    
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
                //TODO: check exact numbering?
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
                ret = snprintf(output_name, MAX_NAME, "%s.txt", argv[optind-1]);
                if (ret >= MAX_NAME || ret < 0) {
                    printf("average_TECprobe_matrices: error - error when storing output name. aborting...\n");
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

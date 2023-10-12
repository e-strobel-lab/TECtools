//
//  merge_values_files.c
//  
//
//  Created by Eric Strobel on 7/26/22.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"

#include "merge_values_files.h"

/* merge_values_files: merge multiple values files into a single file */
//TODO: change vals_cnt variable to something more representative of what it is
int merge_values_files(values_input * vals_ipt, int vals_cnt, char * merged_out_nm, int nonstandard)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int v = 0;  //index of values input file
    int i = 0;  //general purpose index variable
    int j = 0;  //general purpose index variable
    int k = 0;  //general purpose index variable
    
    int got_line[MAX_VALS] = {0};                    //array to flag the success of get_line for each values file
    int proceed = 1;                                 //flag to continue processing
    
    int line_cnt = 0;                                //tracks number of lines in values files
    char line[MAX_VALS][MAX_LINE] = {{0}};           //arrays to store lines from each values file
    char col_nm[XPCTD_FIELDS][MAX_COL_NM+1] = {{0}}; //array to store input column names for error-checking
    char *p_id[MAX_VALS] = {NULL};                   //pointers to variant id
    char *p_vals[MAX_VALS] = {NULL};                 //pointers to variant values
    int field_count = 0;                             //variable for counting number of fields in values file line
    int found_term = 0;                              //flag that terminating null was found
    
    /* generate merged output file */
    FILE * vals_merged = NULL; //pointer for merged output file
    if ((vals_merged = fopen(merged_out_nm, "w")) == NULL) {
        printf("merge_values_files: error - could not open merged output file. Aborting program...\n");
        abort();
    }
    fprintf(vals_merged, "variant_%s", TECdsply_clmn_hdrs[TDSPLY_VID_HDR]); //print id column header
    
    for (line_cnt = 0; proceed; line_cnt++) {   //iterate through each line of the input files
        for (v = 0; v < vals_cnt; v++) {        //process each input file
            
            /* get line from current values file, test that success status is the same
             for all input files. if success status is not the same for all files, the
             input files are not the same length and cannot be processed together. */
            got_line[v] = (get_line(line[v], vals_ipt[v].fp)) ? 1 : 0;  //if successful, got_line=1, else got_line=0
            if (got_line[v] != got_line[0]) { //test if got_line success status equals that of the first input
                printf("merge_values_files: error - input values files are not the same length. aborting...\n");
                abort();
            }
            
            if (got_line[v]) {

                //in each iteration, the number of tabs that precede a non-null character
                //are counted to determine the number of fields in the data line.

                if (line[v][0] == '\t') { //if id field is missing entry, abort
                    printf("merge_values_files: error - values file lines cannot begin with a tab character. aborting...\n");
                    abort();
                    
                } else if (line[v][0]) {  //line contains a at least one field
                    field_count = 1;      //initialize field count to 1
                    
                } else { //it should not be possible to have an empty line here
                    printf("merge_values_files: error - empty line, this should not be possible. aborting...\n");
                    abort();
                }
                
                //split line into id and values strings
                for (i = 0; line[v][i] != '\t' && line[v][i] && i < MAX_LINE; i++) {;} //iterate to first tab
                
                if (line[v][i] == '\t') {       //if the tab separating the id field from the data fields was found
                    line[v][i] = '\0';          //set first tab to null char to split input line
                    p_id[v] = &line[v][0];      //set pointer to variant id string
                    p_vals[v] = &line[v][i+1];  //set pointer to values string
                    
                } else { //unrecognized format error
                    printf("%s\n", line[v]);
                    printf("merge_values_files: error - unrecognized line format. aborting...\n");
                    abort();
                }
                
                
                if (line_cnt == 0) { //reading column header line (first line of each file)
                    
                    if (!nonstandard) { //input data uses standard TECdisplay format, check all column headers
                        
                        if (v == 0) { //reading column header line of first values file
                            
                            /* when reading the column headers of the first values file, parse
                             the headers and validate against the expected headers (variant_id,
                             bound, unbound, fracBound). all subsequent header lines can then
                             be compared diectly against those from the first values file*/
                                                    
                            i = 0;                        //initialize i (column name index) to zero
                            strcpy(col_nm[i++], p_id[0]); //copy the variant id column name to col_nm[0], increment i
                            
                            //parse values headers and store in col_nm array
                            for (j = 0; p_vals[0][j] && i < XPCTD_FIELDS; i++, field_count++) {
                                
                                //copy column header into col_nm array
                                for (k = 0; p_vals[0][j] != '\t' && p_vals[0][j] && k < (MAX_COL_NM); j++, k++) {
                                    col_nm[i][k] = p_vals[0][j];
                                }
                                col_nm[i][k] = '\0';
                                                            
                                //test that loop exited on a tab or null character.
                                //if test fails, column name is too long.
                                if (p_vals[0][j] == '\t') {        //ended on tab, more headers
                                    j++;                           //increment j
                                    
                                } else if (p_vals[0][j] != '\0') { //if loop did not exit on tab or null, name is too long
                                    printf("merge_values_files: error - unexpected long column name in values file. aborting...\n");
                                    abort();
                                }
                            }
                            
                            //test that loop exited on a null character and that
                            //the expected number of value fields were identified
                            if (p_vals[0][j] || field_count != XPCTD_FIELDS) {
                                printf("merge_values_files: error - unexpected headers for values files. aborting...\n");
                                abort();
                            }
                            
                            //test that column headers match expected strings
                            //strstr is used for variant id header because header can be variable but always
                            //contains the substring "id". other headers are not variable, so test is for
                            //an exact match
                            //TODO: make id test look only at last chars of col name?
                            if (strstr(col_nm[TDSPLY_VID_HDR], TECdsply_clmn_hdrs[TDSPLY_VID_HDR]) == NULL ||
                                strcmp(col_nm[TDSPLY_BND_HDR], TECdsply_clmn_hdrs[TDSPLY_BND_HDR])         ||
                                strcmp(col_nm[TDSPLY_UNB_HDR], TECdsply_clmn_hdrs[TDSPLY_UNB_HDR])         ||
                                strcmp(col_nm[TDSPLY_FRC_HDR], TECdsply_clmn_hdrs[TDSPLY_FRC_HDR])) {
                                printf("merge_values_files: error - unexpected headers for values files. aborting...\n");
                                abort();
                            }

                        /* compare all other values file headers to the headers from
                         the first values file, which has already been validated */
                        //TODO: make id test look only at last chars of col name?
                        } else if (strstr(p_id[v], TECdsply_clmn_hdrs[TDSPLY_VID_HDR]) == NULL ||
                                   strcmp(p_vals[v], p_vals[0])) {
                            printf("merge_values_files: error - supplied values files contain discordant headers. aborting...\n");
                            abort();
                        }
                        
                        /* print column headers for the current values file to the merged output file*/
                        fprintf(vals_merged, "\t%s_%s\t%s_%s\t%s_%s",
                                vals_ipt[v].nm, TECdsply_clmn_hdrs[TDSPLY_BND_HDR],
                                vals_ipt[v].nm, TECdsply_clmn_hdrs[TDSPLY_UNB_HDR],
                                vals_ipt[v].nm, TECdsply_clmn_hdrs[TDSPLY_FRC_HDR]);
                        
                    } else { //input data files use non-standard format, only check id column header
                        //TODO: make id test look only at last chars of col name?
                        if (strstr(p_id[v], TECdsply_clmn_hdrs[TDSPLY_VID_HDR]) == NULL) {
                            printf("merge_values_files: error - nonstandard data files must contain the string '%s' in the header of the first column. aborting...\n", TECdsply_clmn_hdrs[TDSPLY_VID_HDR]);
                            abort();
                        }
                        
                        fprintf(vals_merged, "\t%s", p_vals[v]); //print non-standard headers
                    }
                    
                } else { //reading data line
                    
                    /* data lines are validated by comparing the variant id from each input file
                     to the variant id of the first input file. all variant ids must be identical
                     in order for the values files to be merged. */
                    if (strcmp(p_id[v], p_id[0])) {
                        printf("%s\n%s\n", p_id[v], p_id[0]);
                        printf("merge_values_files: error - variant ids are not aligned. aborting...\n");
                        abort();
                        
                    } else if (!nonstandard) { //input data is standard TECdisplay format, check formatting
                        
                        /* read p_vals string to confirm expected number of fields (count started
                         above) and to check that all characters are of expected types*/
                        for (i = 0, found_term = 0; !found_term && i < MAX_LINE; i++) {
                            
                            if (p_vals[v][i] == '\t' || !p_vals[v][i]) { //if a tab or term null was reached
                                
                                if (isdigit(p_vals[v][i-1])) {            //if the preceding character was a digit
                                    field_count++;                       //increment field count
                                }
                                
                                if (!p_vals[v][i]) { //found terminating null
                                    found_term = 1; //set found terminating null flag
                                }
                                
                            } else if (!isdigit(p_vals[v][i]) && //char is not a digit
                                       p_vals[v][i] != '.'    && //or a '.'
                                       p_vals[v][i] != '-') {    //or a '-'
                                
                                printf("merge_values_files: error - unexpected character %c (ASCII: %d) in values string. aborting...\n", p_vals[v][i], p_vals[v][i]);
                                abort();
                            }
                        }
                        
                        if (!found_term) { //if the terminating null was not found, throw error and abort
                            printf("merge_values_file: unexpected long data values line. aborting...\n");
                            abort();
                        }
                        
                        if (field_count != XPCTD_FIELDS) { //check if correct number of fields was read
                            printf("merge_values_files: error - unexpected number of fields in data line. aborting...\n");
                            abort();
                        }
                        
                    }
                    
                    //all data checks passed, or reading non-standard file format
                    if (v == 0) {                            //if reading line from first values file
                        fprintf(vals_merged, "%s", p_id[0]); //print variant id to merged output file
                    }
                    fprintf(vals_merged, "\t%s", p_vals[v]); //print values string to merged output file
                }
            }
        }
        
        if (!got_line[0]) { //reached last line. to reach this part of the code, the end
            proceed = 0;    //of every values file must have been reached simultaneously
        } else {
            fprintf(vals_merged, "\n"); //printed data line, so print newline to start next line
        }
        
    }
    
    /* close merged values file */
    if (fclose(vals_merged) == EOF) {
        printf("merge_values_files: error - error occurred when closing merged input file. Aborting program...\n");
        abort();
    }
    
    return 1;
}

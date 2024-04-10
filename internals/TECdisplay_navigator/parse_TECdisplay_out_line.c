//
//  parse_TECdisplay_out_line.c
//  
//
//  Created by Eric Strobel on 4/4/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../global/global_defs.h"

#include "../TECdisplay_mapper/TECdisplay_output_column_headers.h"

#include "./TECdisplay_navigator_defs.h"
#include "./TECdisplay_navigator_structs.h"


#include "parse_TECdisplay_out_line.h"

void parse_TECdisplay_out_line(char * line, char ** p_id, char ** p_vals, int * bnd, int * unb, double * frc, int mode, int nonstandard)
{
    extern const char TECdsply_clmn_hdrs[4][32]; //column headers from TECdisplay_output_column_headers.c
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    int k = 0; //general purpose index
    
    int ret = 0; //variable for storing snprintf return value
    
    int field_count = 0; //tracks number of fields in line
    int found_term = 0;  //flag that terminating null was found
    
    char tmp_vals[MAX_LINE+1] = {0}; //storage for parsed values string
    
    char col_nm[XPCTD_FIELDS][MAX_COL_NM+1] = {{0}}; //array to store input column names for error-checking
    char *crrnt_val = NULL;                          //pointer to current value field string
    
    if (line[0] == '\t') { //if id field is missing entry, abort
        printf("parse_TECdisplay_out_line: error - values file lines cannot begin with a tab character. aborting...\n");
        abort();
        
    } else if (line[0]) {  //line contains a at least one field
        field_count = 1;   //initialize field count to 1
        
    } else { //it should not be possible to have an empty line here
        printf("parse_TECdisplay_out_line: error - empty line, this should not be possible. aborting...\n");
        abort();
    }
    
    //split line into id and values strings
    for (i = 0; line[i] != '\t' && line[i] && i < MAX_LINE; i++) {;} //iterate to first tab
    
    if (line[i] == '\t') {       //if the tab separating the id field from the data fields was found
        line[i] = '\0';          //set first tab to null char to split input line
        (*p_id) = &line[0];      //set pointer to variant id string
        (*p_vals) = &line[i+1];  //set pointer to values string
        
        ret = snprintf(tmp_vals, MAX_LINE, "%s", *p_vals); //store temporary copy of values string for parsing
        if (ret >= MAX_LINE || ret < 0) {
            printf("parse_TECdisplay_out_line: error - error when storing temporary copy of values string. aborting...");
            abort();
        }
                
    } else { //unrecognized format error
        printf("%s\n", line);
        printf("parse_TECdisplay_out_line: error - unrecognized line format. aborting...\n");
        abort();
    }
    
    if (mode == TDSPLY_HDR_LINE) { //reading column header line
        
        if (!nonstandard) { //input data uses standard TECdisplay format, check all column headers
            
            /* read tmp_vals string to confirm expected number of fields (count started
             above) and check the identity of the column headers */
            
            //parse values headers and store in col_nm array
            i = 0;                      //initialize i to 0
            strcpy(col_nm[i++], *p_id); //store id in index 0 of col_nm array
            for (j = 0; tmp_vals[j] && field_count < XPCTD_FIELDS; i++, field_count++) {
                
                //copy column header into col_nm array
                for (k = 0; tmp_vals[j] != '\t' && tmp_vals[j] && k < MAX_COL_NM; j++, k++) {
                    col_nm[i][k] = tmp_vals[j];
                }
                col_nm[i][k] = '\0'; //terminate column name string
                                           
                //test that loop exited on a tab or null character.
                //if test fails, column name is too long.
                if (tmp_vals[j] == '\t') { //ended on tab, more headers
                    j++;                   //increment j
                    
                } else if (tmp_vals[j] != '\0') { //if loop did not exit on tab or null, name is too long
                    printf("parse_TECdisplay_out_line: error - unexpected long column name in values file. aborting...\n");
                    abort();
                }
            }
            
            //test that loop exited on a null character and that
            //the expected number of value fields were identified
            if (tmp_vals[j] || field_count != XPCTD_FIELDS) {
                printf("parse_TECdisplay_out_line: error - unexpected headers for values files. aborting...\n");
                abort();
            }
            
            //test that column headers match expected strings
            check_header_string(col_nm[TDSPLY_VID_HDR], TECdsply_clmn_hdrs[TDSPLY_VID_HDR]);
            check_header_string(col_nm[TDSPLY_BND_HDR], TECdsply_clmn_hdrs[TDSPLY_BND_HDR]);
            check_header_string(col_nm[TDSPLY_UNB_HDR], TECdsply_clmn_hdrs[TDSPLY_UNB_HDR]);
            check_header_string(col_nm[TDSPLY_FRC_HDR], TECdsply_clmn_hdrs[TDSPLY_FRC_HDR]);
            
        } else { //input data files use non-standard format, only check id column header
            
            if (strstr(*p_id, TECdsply_clmn_hdrs[TDSPLY_VID_HDR]) == NULL) {
                printf("parse_TECdisplay_out_line: error - nonstandard data files must contain the string '%s' in the header of the first column. aborting...\n", TECdsply_clmn_hdrs[TDSPLY_VID_HDR]);
                abort();
            }
        }
        
    } else if (mode == TDSPLY_DATA_LINE) { //reading data line
        
        if (!nonstandard) { //input data is standard TECdisplay format, check formatting
            
            /* read tmp_vals string to confirm expected number of fields (count started
             above) and to check that all characters are of expected types */
            
            //found_term flag is used to include processing of the last value field in the loop
            for (i = 0, crrnt_val = &tmp_vals[0], found_term = 0; !found_term && i < MAX_LINE; i++) {
                
                if (tmp_vals[i] == '\t' || !tmp_vals[i]) { //if a tab or term null was reached
                    
                    if (!tmp_vals[i]) {  //found terminating null
                        found_term = 1;  //set found terminating null flag
                    }
                    
                    tmp_vals[i] = '\0';  //set the delimiter to a term null
                    
                    if (field_count == TDSPLY_BND_HDR) {         //if reading bound reads field
                        *bnd = atoi(crrnt_val);                  //store bound read count
                    
                    } else if (field_count == TDSPLY_UNB_HDR) {  //if reading unbound reads field
                        *unb = atoi(crrnt_val);                  //store unbound read count
                        
                    } else if (field_count == TDSPLY_FRC_HDR) {  //if reading frac bound reads field
                        *frc = strtof(crrnt_val, NULL);          //store frac bound value
                    }
                    
                    field_count++;  //increment field count
                    
                    crrnt_val = &tmp_vals[i+1];  //set pointer to start of next value field
                    
                } else { //verify that character is allowable in field string
                    
                    //bound and unbound read counts are integers and allow only digits
                    if (field_count == TDSPLY_BND_HDR || field_count == TDSPLY_UNB_HDR) {
                        if (!isdigit(tmp_vals[i])) {
                            printf("parse_TECdisplay_out_line: error - unexpected non-digit character %c (ASCII: %d) in integer values string. aborting...\n", tmp_vals[i], tmp_vals[i]);
                            abort();
                        }
                        
                    //fraction bound allows float chars
                    } else if (field_count == TDSPLY_FRC_HDR) {
                        if (!isdigit(tmp_vals[i]) && //char is not a digit
                                   tmp_vals[i] != '.'    && //or a '.'
                                   tmp_vals[i] != '-'    && //or a '-'
                                   tmp_vals[i] != 'e'    && //or an 'e'
                                   tmp_vals[i] != 'n'    && //or an 'n'
                                   tmp_vals[i] != 'a'    ){ //or an 'a'
                            printf("parse_TECdisplay_out_line: error - unexpected character %c (ASCII: %d) in float values string. aborting...\n", tmp_vals[i], tmp_vals[i]);
                            abort();
                        }
                    }
                }
            }
            
            if (!found_term) { //if the terminating null was not found, throw error and abort
                printf("parse_TECdisplay_out_line: unexpected long data values line. aborting...\n");
                abort();
            }
            
            if (field_count != XPCTD_FIELDS) { //check if correct number of fields was read
                printf("parse_TECdisplay_out_line: error - unexpected number of fields in data line. aborting...\n");
                abort();
            }
        }
        
    } else {
        printf("parse_TECdisplay_out_line: error - unrecognized mode. aborting...\n");
        abort();
    }

}

void check_header_string(char * haystack, const char * needle)
{
    int len = strlen(haystack); //string length/terminating null index
    
    char *p_test = NULL; //pointer to string to test
    
    int i = 0;           //general purpose index
    int fnd_uscore = 0;  //flag that underscore was found
    
    i = len;                         //set i to terminating null index
    while (i >= 0 && !fnd_uscore) {  //until the start of the string or an underscore is reached
        if (haystack[i] == '_') {    //if an underscore is found
            fnd_uscore = 1;          //set flag that the underscore was found
        } else {                     //otherwise
            i--;                     //decrement i
        }
    }
    
    if (fnd_uscore) {             //if an underscore was found
        p_test = &haystack[i+1];  //set test pointer to the char after the underscore
        
    } else {                      //if no underscore was found
        p_test = &haystack[0];    //set test pointer to start of haystack string
    }
    
    if (strcmp(p_test, needle)) {  //compare test pointer to needle
        printf("parse_TECdisplay_out_line: error - did not find '%s' column header in expected location. aborting...\n", needle);
        abort();
    }
}

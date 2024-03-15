//
//  count_delims_2_col.c
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#include "../global/global_defs.h"
#include "./assemble_TECprobeLMP_data_defs.h"
#include "./assemble_TECprobeLMP_data_structs.h"

#include "count_delims_2_col.h"

/* count_delims_2_col: count the number of delimiters observed before
 the target column header is reached */
int count_delims_2_col(FILE * ipt, mode_parameters * mode_params, int nrchd_len, int * delims2col)
{
    int i = 0;                            //general purpose index
    int j = 0;                            //general purpose index
    
    char c = '\0';                        //current character of the line
     
    int found_end = 0;                    //flag that end of line was found
    
    char crnt_col_nm[MAX_NAME+1] = {0}; //array to store column name TODO: allocate dynamically?
    
    char * trgt_hdr_strt = NULL;          //pointer to start of target header string
    char * trgt_len_strt = NULL;          //pointer to start of length substring in target header string (REACTIVITY mode)
    
    for (i = 0; !found_end; i++) {        //until the end of the line is found
        
        crnt_col_nm[0] = '\0';            //set crnt_col_nm to null char
        trgt_hdr_strt = NULL;             //set trgt_hdr_strt pointer to NULL
        trgt_len_strt = NULL;             //set trgt_len_strt pointer to NULL
         
        //get characters until a delimiter or the end of the line is reached
        //and store the string as the column name
        for (j = 0; (c = fgetc(ipt)) != mode_params->dlm && j < MAX_NAME && c != '\n' && c != EOF; j++) {
            crnt_col_nm[j] = c;
        }
        crnt_col_nm[j] = '\0';
        
        if (j == MAX_NAME) { //check that current column name does not exceed maximum length
            printf("count_delims_2_col: error - column name starting with %s is too long. aborting...\n", crnt_col_nm);
            abort();
        }
        
        if (c == '\n' || c == EOF) { //if loop ended on a newline or EOF
            found_end = 1;           //set flag that the end of the line was found
        }
        
        if ((trgt_hdr_strt = strstr(crnt_col_nm, mode_params->hdr)) != NULL) { //search column header for target header substring
            
            if (mode_params->mod == REACTIVITY) {                         //if running REACTIVITY mode
                trgt_len_strt = &trgt_hdr_strt[strlen(mode_params->hdr)]; //set pointer to the start of the length substring in target header
            }
            
            //in REACTIVITY mode, check that the current column is a match to the current enriched length
            //in LENGTH_DISTRIBUTION mode, the match is already established so always proceed
            if ((mode_params->mod == REACTIVITY && atoi(trgt_len_strt) == nrchd_len) ||
                (mode_params->mod == LEN_DIST)) {
                
                *delims2col = i; //set the number of delimiters that precede the target column
                
                if (c != '\n' && c != EOF) {                           //if the current line character is not a newline or EOF
                    while ((c = fgetc(ipt)) != '\n' && c != EOF) { ; } //iterate to the end of the line
                }
                
                return 1; //return success
            }
        }
    }
    
    return 0; //return failure
}

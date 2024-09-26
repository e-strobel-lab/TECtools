//
//  get_value.c
//  
//
//  Created by Eric Strobel on 2/7/24.
//

#include "../global/global_defs.h"
#include "./assemble_TECprobeLM_data_defs.h"
#include "./assemble_TECprobeLM_data_structs.h"

#include "get_value.h"

/* get_value: get the value of the target column in the current line */
int get_value(FILE * ipt,  mode_parameters * mode_params, int delims2col, int delims2seq, double * val, char ** val_str, char * seq)
{
    int i = 0;                         //general purpose index
    int j = 0;                         //general purpose index
    int k = 0;                         //general purpose index
    
    char c = '\0';                     //current character of the line
    
    char tmp_val[MAX_NAME+1] = {0};    //storage for column value string
    int delim_cnt = 0;                 //delimiter counter
    int fnd_trgt_col = 0;              //flag that the target column was found
    
    int set_val = 0;                   //flag that target value was set
    int set_seq = 0;                   //flag that sequence was set
    
    for (i = 0, k = 0, c = '\0', fnd_trgt_col = 0; !fnd_trgt_col && c != '\n' && c != EOF; i++) {
        
        c = fgetc(ipt);       //get the next character
        
        if (k < MAX_NAME) { //check that tmp_val array bounds was not exceeded
            tmp_val[k++] = c; //store the current character in tmp_val
        } else {
            printf("get_value: error - column value starting with %s is too long. aborting\n", tmp_val);
            abort();
        }
        
        if (c == mode_params->dlm) { //if current character is a delimiter
            delim_cnt++;             //increment the delimiter counter
            
            tmp_val[k-1] = '\0';     //terminate the value string
            k = 0;                   //reset k to 0
            
            //in LENGTH_DISTRIBUTION mode, if offset has not been set yet and
            //the first column value (transcript length) has been read, set the
            //offset to the current column value (which is the minimum transcript
            //length)
            if (mode_params->mod == LEN_DIST && !mode_params->offset && delim_cnt == 1) {
                mode_params->offset = atoi(tmp_val);
            }
            
            if (delim_cnt == delims2col) { //if the target column has been reached
                for (j = 0; (c = fgetc(ipt)) != mode_params->dlm && j < MAX_NAME && c != '\n' && c != EOF; j++) {
                    tmp_val[j] = c;
                }
                tmp_val[j] = '\0';
                
                if (j < MAX_NAME) {    //check that tmp_val array bounds was not exceeded
                    *val = strtod(tmp_val, NULL); //store the value of the current column
                    if (mode_params->mod == REACTIVITY) {
                        if (((*val_str) = malloc((strlen(tmp_val)+1) * sizeof(**val_str))) == NULL) {
                            printf("get_value: error - memory allocation for matrix field value failed. aborting...\n");
                            abort();
                        } else {
                            strcpy(*val_str, tmp_val); //copy tmp_val to mtrx-vals array
                        }
                    }
                    set_val = 1; //change flag that value was set to true
                    
                } else {
                    printf("get_value: error - column value starting with %s is too long. aborting\n", tmp_val);
                    abort();
                }
                
            } else if (mode_params->mod == REACTIVITY && delim_cnt == delims2seq) { //store target RNA sequence
                
                if (strlen(tmp_val) == 1) {           //check that sequence is a single character
                    if (!(*seq)) {                    //if sequence was not yet set
                        *seq = tmp_val[0];            //set sequence
                    } else if (tmp_val[0] != *seq) {  //otherwise check that character matches previous sequence
                        printf("get_value: error - discordant sequences detected in input files. aborting...\n");
                        abort();
                    }
                } else {
                    printf("get_value: error - sequence string is more than one character. aborting...\n");
                    abort();
                }
                set_seq = 1; //change flag that sequence was set to true
            }
            
            if ((mode_params->mod == REACTIVITY && (set_val && set_seq)) ||
                (mode_params->mod == LEN_DIST && set_val)) {
                
                while (c != '\n' && c != EOF) { //iterate to the end of the line
                    c = fgetc(ipt);
                }
                
                return c; //return the terminating character of the line
            }
        }
    }
    
    if (i == 1 && c == EOF) { //if the first character of the line was EOF
        return c;             //return c (== EOF)
    } else {
        return 0;             //return failue
    }
}

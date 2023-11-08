//
//  get_constraint_metadata.c
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#include <stdio.h>
#include <ctype.h>

#include "../global/global_defs.h"

#include "../utils/io_management.h"

#include "../seq_utils/basemap.h"

#include "../TECdisplay_navigator/parse_reference.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

#include "get_constraint_metadata.h"

void get_constraint_metadata(char * ipt_fn, constraint_metadata * cons_meta, char type)
{
    get_file(&(cons_meta->fp), ipt_fn);     //set file pointer to constraints file
    strcpy(cons_meta->fn, ipt_fn);          //store file name
    get_sample_name(ipt_fn, cons_meta->sn); //get sample name from constraints file name
    cons_meta->typ = type;                  //set constraint type
    
    FILE * tmp_ifp; //temporary file pointer for getting constraint names
    
    //int i = 0; //index for constraint groups
    int j = 0; //index for input line
    int k = 0; //index for code and tmp arrays
    int n = 0; //index for checking duplicate constraint names
    
    char line[MAX_LINE+1] = {0};  //array to store input line
    char code[MAX_CODE+1] = {0};  //array to store constraint code
    
    int new_constraint = 1;       //flag to indicate new constraint
    
    get_file(&(tmp_ifp), ipt_fn); //set temp file pointer to constraints file
        
    parse_reference(tmp_ifp, &cons_meta->bmap, &cons_meta->wt, &cons_meta->cnstnt_indels);
    
    while (get_line(line, tmp_ifp)) {
                
        if (new_constraint) {
            /***** parsing the first line of a new constraint, store constraint name *****/
            
            //check that number of constraints doesn't exceed MAX_CONSTRAINTS
            if (cons_meta->c_cnt >= MAX_CONSTRAINTS) {
                printf("get_constraint_metadata: error - number of supplied constraints exceeds the expected maximum (%d). aborting...\n", MAX_CONSTRAINTS);
                abort();
            }
            
            //copy name to constraints struct, replace any spaces with underscores
            //newlines are trimmed by get_line, so only internal spaces should remain
            //replacing spaces with underscores is necessary for datagraph compatibility
            for (j = 0; line[j] && j < MAX_LINE; j++) {
                cons_meta->cn[cons_meta->c_cnt][j] = (isspace(line[j])) ? '_' : line[j];
            }
            cons_meta->cn[cons_meta->c_cnt][j] = '\0';
            
            //check for duplicate constraint name
            for (n = 0; n < cons_meta->c_cnt; n++) {
                if (!strcmp(cons_meta->cn[n], cons_meta->cn[cons_meta->c_cnt])) {
                    printf("parse_constraints: error - duplicate entry for constraint %s. aborting...\n", cons_meta->cn[cons_meta->c_cnt]);
                    abort();
                }
            }
            
            cons_meta->c_cnt++;
            new_constraint = 0; //turn off new constraint flag
            
        } else {
            /********* parsing body of a constraint *********/
            
            code[0] = '\0'; //zero code arrays

            /* get constraint code
             base = constrained base
             pair = constrained pair
             # = end of constraint*/
            for (j = 0, k = 0; !isspace(line[j]) && line[j] && k < MAX_CODE; j++, k++) {
                code[k] = line[j];
            }
            code[k] = '\0';
            
            if (!isspace(line[j]) && code[0] != '#') {
                printf("parse_constraints: error - unexpected format for constraint line.\n%s\naborting...\n", line);
                abort();
            }
                        
            /********* process line based on constraint code *********/
            if (code[0] == '#') { //reached end of constraint

                new_constraint = 1; //turn on new constraint flag
                
            } else if (!strcmp(code, "base")) { //line specifies constrained base
                ;
            } else if (!strcmp(code, "pair")) { //line specifies constrained pair
                ;
            } else { //unrecognized constraint code
                printf("parse_constraints: unrecognized constraint code \"%s\". aborting...", code);
                abort();
            }
        }
    }
    
    /* close output file */
    if (fclose(tmp_ifp) == EOF) {
        printf("filter_values: error - error occurred when closing output file %s.txt. Aborting program...\n", ipt_fn);
        abort();
    }
}

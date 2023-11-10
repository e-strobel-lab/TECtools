//
//  calc_output_files.c
//  
//
//  Created by Eric Strobel on 11/7/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../utils/io_management.h"

#include "./TECdisplay_Hnav_defs.h"
#include "./TECdisplay_Hnav_structs.h"
#include "./TECdisplay_Hnav_global_vars.h"

#include "calc_output_files.h"

/* calc_output_files: calculate the expected number of output files and
 confirm that TECdisplay_navigator analysis should proceed */
int calc_output_files(int layr_cnt, constraint_metadata * cons_meta)
{
    extern output_file_names out_fns[MAX_LAYERS]; //storage for output file names of each layer
    
    int i = 0; //general purpose index
    int j = 0; //general purpose index
    
    int f_cnt = 1;     //file count
    int tot_files = 0; //total files that will be generated
    
    char usr_resp[4] = {0};       //storage for user response
    char discard[MAX_LINE] = {0}; //array for flushing stdin
    int resp_provided = 0;        //flag that valid reponse was provided
    
    printf("the provided constraints will generate the following files:\n\n");
    
    for (i = 0; i < layr_cnt; i++) {     //for each layer
        if (cons_meta[i].typ == 'c') {   //if constraint matches will be included
            f_cnt *= cons_meta[i].c_cnt; //multiply the file count by the number of constraints in the constraint file
        }                                //if constraint matches are excluded, f_cnt will be the same as the previous layer
        
        out_fns[i].f_cnt = f_cnt; //store the file count for the current layer
        tot_files += f_cnt;       //increment the total file count by the file count of the current layer
        
        printf("layer %d: %d files\n", i+1, f_cnt); //print the file count for the current layer
        
        //allocate pointers for file names
        if ((out_fns[i].fn = calloc(f_cnt, sizeof(*(out_fns[i].fn)))) == NULL) {
            printf("calc_output_files: error - memory allocation output file names failed. aborting...\n");
            abort();
        }
    }
    
    printf("in total, %d files will be generated.\n\n", tot_files); //print the total file count
    
    printf("proceed? (yes/no)\n");
    
    while (!resp_provided) {
        scanf("%3s", usr_resp);
        
        if (!strcmp(usr_resp, "yes")) {
            printf("proceeding with analysis\n");
            resp_provided = 1;
            
        } else if (!strcmp(usr_resp, "no")) {
            printf("aborting analysis\n");
            abort();
            
        } else {
            printf("invalid response. please enter \"yes\" or \"no\"\n");
            get_line(discard, stdin);
            
        }
    }
    
    return tot_files;
}
